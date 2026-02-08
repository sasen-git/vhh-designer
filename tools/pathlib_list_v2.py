#!/usr/bin/env python3
from pathlib import Path
from datetime import datetime

# Optional: folders to skip
EXCLUDE_DIRS = {
    ".git", "__pycache__", ".venv", "venv", "node_modules",
    ".idea", ".vscode", "dist", "build"
}

# Where to write the output file (no new folder is created under the scanned root)
OUT_DIR = Path("/home/sasenefrem/KA-Search/active/utilities/file_list")

def main() -> None:
    root = Path.cwd()

    # Ensure output directory exists
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Date + time in filename (e.g., file_list_2026-01-02_01-23-45.txt)
    dt_str = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    out_file = OUT_DIR / f"file_list_{dt_str}.txt"
    out_file_resolved = out_file.resolve()

    lines = []
    lines.append("FILE LIST")
    lines.append(f"Root: {root.resolve()}")
    lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")

    out_dir_resolved = OUT_DIR.resolve()

    for path in sorted(root.rglob("*")):
        # Skip excluded directories (and anything inside them)
        if any(part in EXCLUDE_DIRS for part in path.parts):
            continue

        # If the output directory is inside the scanned root, skip it and its contents
        try:
            path.resolve().relative_to(out_dir_resolved)
            continue
        except ValueError:
            pass

        # Skip the output file itself (in case it lands inside the scanned root)
        if path.is_file() and path.resolve() == out_file_resolved:
            continue

        rel = path.relative_to(root)

        if path.is_dir():
            lines.append(f"[D] {rel.as_posix()}/")
        elif path.is_file():
            lines.append(f"[F] {rel.as_posix()}")

    out_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote: {out_file}")

if __name__ == "__main__":
    main()
