from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class DataSourceConfig:
    """Archive access settings for one telescope or survey family."""

    status: str
    archive: str
    notes: tuple[str, ...] = ()
    starting_collections: tuple[str, ...] = ()


@dataclass(frozen=True)
class CandidateGenerationConfig:
    """First-pass residual detection settings."""

    threshold_sigma: float = 5.0
    min_epochs: int = 3
    spatial_radius_arcsec: float = 0.08
    temporal_window_days: float = 120.0


@dataclass(frozen=True)
class ModelConfig:
    """High-level model interface."""

    name: str
    classes: tuple[str, ...]
    image_shape: tuple[int, int, int, int] = (12, 5, 31, 31)


@dataclass(frozen=True)
class ArchitectureConfig:
    """Typed subset of the YAML architecture configuration."""

    project_name: str
    science_goal: str
    roman: DataSourceConfig
    hubble: DataSourceConfig
    candidate_generation: CandidateGenerationConfig
    model: ModelConfig
    raw: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "ArchitectureConfig":
        raw = _load_yaml(Path(path))

        project = raw["project"]
        data_sources = raw["data_sources"]
        candidate = raw["candidate_generation"]
        detection = candidate.get("detection", {})
        grouping = candidate.get("grouping", {})
        model = raw["model"]
        image_shape = tuple(model["inputs"]["image_sequence"]["shape"])

        return cls(
            project_name=project["name"],
            science_goal=project["science_goal"],
            roman=_source_from_dict(data_sources["roman"]),
            hubble=_source_from_dict(data_sources["hubble"]),
            candidate_generation=CandidateGenerationConfig(
                threshold_sigma=float(detection.get("threshold_sigma", 5.0)),
                min_epochs=int(detection.get("min_epochs", 3)),
                spatial_radius_arcsec=float(grouping.get("spatial_radius_arcsec", 0.08)),
                temporal_window_days=float(grouping.get("temporal_window_days", 120.0)),
            ),
            model=ModelConfig(
                name=model["name"],
                classes=tuple(model["outputs"]["classes"]),
                image_shape=image_shape,  # type: ignore[arg-type]
            ),
            raw=raw,
        )


def _source_from_dict(raw: dict[str, Any]) -> DataSourceConfig:
    return DataSourceConfig(
        status=raw["status"],
        archive=raw["archive"],
        notes=tuple(raw.get("notes", ())),
        starting_collections=tuple(raw.get("starting_collections", ())),
    )


def _load_yaml(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    try:
        import yaml  # type: ignore[import-not-found]
    except ImportError:
        return _load_minimal_yaml(text)
    return yaml.safe_load(text)


def _load_minimal_yaml(text: str) -> dict[str, Any]:
    """Parse the small YAML subset used by this scaffold.

    PyYAML is the right parser for production. This fallback keeps the example
    usable in a brand-new Python environment where optional dependencies have
    not been installed yet.
    """

    lines = [
        (len(line) - len(line.lstrip(" ")), line.strip())
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    parsed, index = _parse_yaml_block(lines, 0, 0)
    if index != len(lines):
        raise ValueError("Could not parse full YAML configuration")
    if not isinstance(parsed, dict):
        raise ValueError("Top-level YAML configuration must be a mapping")
    return parsed


def _parse_yaml_block(
    lines: list[tuple[int, str]], index: int, indent: int
) -> tuple[Any, int]:
    if index >= len(lines):
        return {}, index

    current_indent, content = lines[index]
    if current_indent < indent:
        return {}, index
    if current_indent != indent:
        raise ValueError(f"Unexpected indentation near: {content}")

    if content.startswith("- "):
        return _parse_yaml_list(lines, index, indent)
    return _parse_yaml_mapping(lines, index, indent)


def _parse_yaml_mapping(
    lines: list[tuple[int, str]], index: int, indent: int
) -> tuple[dict[str, Any], int]:
    result: dict[str, Any] = {}
    while index < len(lines):
        current_indent, content = lines[index]
        if current_indent < indent:
            break
        if current_indent > indent:
            raise ValueError(f"Unexpected nested value near: {content}")
        if content.startswith("- "):
            break
        if ":" not in content:
            raise ValueError(f"Expected key/value pair near: {content}")

        key, raw_value = content.split(":", 1)
        raw_value = raw_value.strip()
        if raw_value:
            result[key] = _parse_yaml_scalar(raw_value)
            index += 1
            continue

        if index + 1 >= len(lines) or lines[index + 1][0] <= indent:
            result[key] = {}
            index += 1
            continue

        result[key], index = _parse_yaml_block(lines, index + 1, lines[index + 1][0])
    return result, index


def _parse_yaml_list(
    lines: list[tuple[int, str]], index: int, indent: int
) -> tuple[list[Any], int]:
    result: list[Any] = []
    while index < len(lines):
        current_indent, content = lines[index]
        if current_indent < indent:
            break
        if current_indent != indent or not content.startswith("- "):
            break

        raw_item = content[2:].strip()
        if raw_item:
            result.append(_parse_yaml_scalar(raw_item))
            index += 1
            continue

        if index + 1 >= len(lines):
            result.append({})
            index += 1
            continue
        item, index = _parse_yaml_block(lines, index + 1, lines[index + 1][0])
        result.append(item)
    return result, index


def _parse_yaml_scalar(value: str) -> Any:
    if value.startswith("[") and value.endswith("]"):
        raw_items = value[1:-1].strip()
        if not raw_items:
            return []
        return [_parse_yaml_scalar(item.strip()) for item in raw_items.split(",")]

    if value in {"true", "True"}:
        return True
    if value in {"false", "False"}:
        return False
    if value in {"null", "None"}:
        return None

    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass

    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    return value
