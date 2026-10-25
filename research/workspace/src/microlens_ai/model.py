from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HybridModelSpec:
    """Framework-neutral model description."""

    image_embedding_dim: int = 256
    light_curve_embedding_dim: int = 256
    physics_feature_dim: int = 16
    fusion_dim: int = 384
    class_count: int = 7
    regression_count: int = 4


def describe_model(spec: HybridModelSpec | None = None) -> dict[str, object]:
    """Return a serializable model plan for dashboards and docs."""

    spec = spec or HybridModelSpec()
    return {
        "name": "hybrid_microlensing_transformer",
        "branches": {
            "image": "2D CNN per epoch plus temporal attention",
            "light_curve": "irregular-time transformer with uncertainty features",
            "physics": "analytic microlensing and variable-template fit features",
        },
        "fusion": "gated attention over image, light-curve, and physics embeddings",
        "outputs": {
            "classes": spec.class_count,
            "regression": spec.regression_count,
            "uncertainty": ["epistemic", "aleatoric"],
        },
    }


class HybridMicrolensingNet:
    """Optional PyTorch implementation shell.

    Importing this module does not require torch. Instantiating the class does.
    """

    def __init__(self, spec: HybridModelSpec | None = None) -> None:
        try:
            import torch
            from torch import nn
        except ImportError as exc:
            raise RuntimeError(
                "Install the optional 'ml' dependencies to instantiate the PyTorch model."
            ) from exc

        self.spec = spec or HybridModelSpec()
        self.torch = torch
        self.nn = nn
        self.image_encoder = nn.Sequential(
            nn.Conv2d(5, 32, kernel_size=3, padding=1),
            nn.GELU(),
            nn.Conv2d(32, 64, kernel_size=3, padding=1),
            nn.GELU(),
            nn.AdaptiveAvgPool2d((1, 1)),
            nn.Flatten(),
            nn.Linear(64, self.spec.image_embedding_dim),
        )
        self.light_curve_encoder = nn.Sequential(
            nn.Linear(8, 128),
            nn.GELU(),
            nn.Linear(128, self.spec.light_curve_embedding_dim),
        )
        fusion_input = (
            self.spec.image_embedding_dim
            + self.spec.light_curve_embedding_dim
            + self.spec.physics_feature_dim
        )
        self.fusion = nn.Sequential(
            nn.Linear(fusion_input, self.spec.fusion_dim),
            nn.GELU(),
            nn.Dropout(0.1),
        )
        self.class_head = nn.Linear(self.spec.fusion_dim, self.spec.class_count)
        self.regression_head = nn.Linear(self.spec.fusion_dim, self.spec.regression_count)

    def parameters(self):
        modules = [
            self.image_encoder,
            self.light_curve_encoder,
            self.fusion,
            self.class_head,
            self.regression_head,
        ]
        for module in modules:
            yield from module.parameters()
