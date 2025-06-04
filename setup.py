import setuptools

with open("README.md", "r", encoding = "utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "repex_topology_parser",
    version = "0.1.0",
    author = "Korey M. Reid",
    author_email = "r.korey@gmail.com",
    description = "Parse gromacs topology for scaling",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/koreyr/repex_topology_parser.git",
    project_urls = {
        "Bug Tracker": "https://github.com/koreyr/repex_topology_parser/issues",
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir = {"": "src"},
    packages = setuptools.find_packages(where="src"),
    python_requires = ">=3.6",
    install_requires=[
        'pandas >=2.2',
        'numpy <2',
    ],
)
