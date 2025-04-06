#!/usr/bin/env python3

import argparse
import logging
import subprocess
import shutil
from pathlib import Path

# Configure logging
LOGGER = logging.getLogger("git_subclone")
handler = logging.StreamHandler()
formatter = logging.Formatter(
    fmt="%(asctime)s %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
handler.setFormatter(formatter)
LOGGER.addHandler(handler)
LOGGER.setLevel(logging.INFO)

# List of repos to clone in the format (repo_url, commit_sha_or_tag, subdirectory_name)
REPO_SERIES = [
    (
        "git@gitlab.internal.sanger.ac.uk:DERMATLAS/fur/fur_germline.git",
        "0.3.1",
        "germline-plot",
    ),
    (
        "git@gitlab.internal.sanger.ac.uk:DERMATLAS/fur/fur2cosmic.git",
        "477d7a31921c47fe9779a53e260376ce4318e949",
        "human-comparison-plot",
    ),
]


def clone_repo(repo_url, ref, subdir_path, force=False):
    """
    Clone a git repository, checkout a specific reference, and remove the .git directory.

    Args:
        repo_url (str): The URL of the git repository to clone
        ref (str): The commit SHA or tag to checkout
        subdir_path (Path): Path to the subdirectory where the repo will be cloned
        force (bool): If True, delete the subdirectory if it already exists

    Returns:
        bool: True if successful, False otherwise
    """
    # Check if the subdirectory already exists
    if subdir_path.exists():
        if force:
            LOGGER.info(f"Removing existing directory: {subdir_path}")
            try:
                shutil.rmtree(subdir_path)
            except Exception as e:
                LOGGER.error(f"Failed to remove directory: {subdir_path}")
                LOGGER.error(f"Error: {e}")
                return False
        else:
            LOGGER.info(
                f"Directory already exists: {subdir_path}. Skipping. Use --force to overwrite."
            )
            return True

    # Clone the repository
    LOGGER.info(f"Cloning {repo_url} into {subdir_path}")
    result = subprocess.run(
        ["git", "clone", repo_url, str(subdir_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        LOGGER.error(f"Failed to clone repository: {repo_url}")
        LOGGER.error(f"Error: {result.stderr}")
        return False

    # Checkout the specified reference
    LOGGER.info(f"Checking out {ref}")
    result = subprocess.run(
        ["git", "checkout", ref],
        cwd=str(subdir_path),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        LOGGER.error(f"Failed to checkout reference: {ref}")
        LOGGER.error(f"Error: {result.stderr}")
        return False

    # Remove the .git directory
    git_dir = subdir_path / ".git"
    if git_dir.exists():
        LOGGER.info(f"Removing .git directory from {subdir_path}")
        shutil.rmtree(git_dir)

    LOGGER.info(f"Successfully cloned {repo_url} at {ref} into {subdir_path}")
    return True


def main():
    """
    Main function to parse arguments and clone repositories.

    Returns:
        int: 0 for success, 1 for failure
    """
    parser = argparse.ArgumentParser(
        description="Clone git repositories into subdirectories and remove .git directories"
    )
    parser.add_argument(
        "--force", action="store_true", help="Force overwrite existing directories"
    )
    parser.add_argument(
        "base_dir",
        nargs="?",
        default=str(Path.cwd()),
        help="Base directory for cloning repositories (default: current working directory)",
    )

    args = parser.parse_args()

    # Convert base_dir to a Path object
    base_dir = Path(args.base_dir)

    # Ensure the base directory exists
    base_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info(f"Using base directory: {base_dir}")

    # Clone each repository
    success = True
    for repo_url, ref, subdir_name in REPO_SERIES:
        subdir_path = base_dir / subdir_name
        if not clone_repo(repo_url, ref, subdir_path, args.force):
            success = False

    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
