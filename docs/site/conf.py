"""Documentation configuration; building this site launches no proof checks."""

project = "Coxeter Polytopes"
author = "CoxeterPolytopes contributors"
copyright = "2026, CoxeterPolytopes contributors"
language = "en"
root_doc = "index"
extensions = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "alabaster"
html_title = "Coxeter Polytopes: papers and verification"
html_baseurl = "https://sashakolpakov.github.io/CoxeterPolytopes/"
html_static_path = []
templates_path = ["_templates"]
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        "sourcelink.html",
    ],
}
html_theme_options = {
    "description": "Six papers, explicit proof inputs, and bounded verification",
    "github_user": "sashakolpakov",
    "github_repo": "CoxeterPolytopes",
    "github_banner": True,
    "github_button": False,
}
html_context = {
    "source_repository": "https://github.com/sashakolpakov/CoxeterPolytopes",
    "source_branch": "main",
    "source_directory": "docs/site",
}

