from ford.fortran_project import Project
from ford.settings import load_toml_settings, load_markdown_settings


def get_project() -> Project:
    directory = "."

    # load settings from fpm.toml
    proj_settings = load_toml_settings(directory)

    # if no fpm.toml, use ford.yml
    if proj_settings is None:
        with open("ford.yml", "r") as f:
            proj_settings, _ = load_markdown_settings(directory, f.read(), f.name)

    return Project(proj_settings)
