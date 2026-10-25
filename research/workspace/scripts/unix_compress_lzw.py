"""Small UNIX compress (.Z) LZW decoder used for ESO archive FITS streams.

The ESO data portal can return classic UNIX compress payloads with a `.fits`
suffix. Python's standard library does not decode this format, and the local
uncompresspy install in this workspace is incomplete. This helper decodes the
payload into bytes so astropy can open the FITS file from memory.
"""

from __future__ import annotations

from pathlib import Path


def decompress_z_bytes(payload: bytes) -> bytes:
    if len(payload) < 3 or payload[0] != 0x1F or payload[1] != 0x9D:
        return payload

    flags = payload[2]
    maxbits = flags & 0x1F
    block_mode = bool(flags & 0x80)
    if maxbits < 9 or maxbits > 16:
        raise ValueError(f"Unsupported .Z maxbits={maxbits}")

    bitpos = 24
    width = 9
    clear_code = 256
    next_code = 257 if block_mode else 256
    table = {i: bytes([i]) for i in range(256)}
    output = bytearray()
    old_code: int | None = None

    def reset() -> None:
        nonlocal width, next_code, table, old_code
        width = 9
        next_code = 257 if block_mode else 256
        table = {i: bytes([i]) for i in range(256)}
        old_code = None

    def read_code(nbits: int) -> int | None:
        nonlocal bitpos
        if bitpos + nbits > len(payload) * 8:
            return None
        value = 0
        for i in range(nbits):
            byte = payload[(bitpos + i) >> 3]
            bit = (byte >> ((bitpos + i) & 7)) & 1
            value |= bit << i
        bitpos += nbits
        return value

    while True:
        code = read_code(width)
        if code is None:
            break

        if block_mode and code == clear_code:
            reset()
            continue

        if old_code is None:
            if code not in table:
                break
            entry = table[code]
            output.extend(entry)
            old_code = code
            continue

        if code in table:
            entry = table[code]
        elif code == next_code:
            previous = table[old_code]
            entry = previous + previous[:1]
        else:
            break

        output.extend(entry)
        if next_code < (1 << maxbits):
            table[next_code] = table[old_code] + entry[:1]
            next_code += 1
            if next_code >= (1 << width) and width < maxbits:
                width += 1
        old_code = code

    return bytes(output)


def decompress_z_file(path: str | Path) -> bytes:
    return decompress_z_bytes(Path(path).read_bytes())

