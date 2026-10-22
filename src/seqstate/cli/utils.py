from pathlib import Path


def create_versioned_path(filepath: Path) -> Path:
    if not filepath.parent.is_dir():
        filepath.parent.mkdir(parents=True, exist_ok=True)

    version = 0
    while (
        path := Path(filepath.parent, f"{filepath.stem}_v{version}{filepath.suffix}")
    ).exists():
        version += 1

    return path
