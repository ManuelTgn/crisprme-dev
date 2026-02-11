from setuptools import setup, find_packages
from setuptools_rust import RustExtension, Binding

# Read README for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [
        line.strip() for line in fh if line.strip() and not line.startswith("#")
    ]

setup(
    name="crisprme2",
    version="0.1.0",
    author="Manuel Tognon",
    author_email="manu.tognon@gmail.com",
    maintainer="Manuel Tognon",
    maintainer_email="manu.tognon@gmail.com",
    description=(
        "CRISPRme2: High-performance and scalable tool for variant- and "
        "haplotype-aware genome-wide off-target assessment in CRISPR-Cas systems"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ManuelTgn/CRISPRme2",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    extras_require={
        "dev": ["pytest", "black", "flake8", "mypy"],
        "docs": ["sphinx", "sphinx-rtd-theme"],
    },
    entry_points={
        "console_scripts": [
            "crisprme2=crisprme2.__main__:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    # ===== Rust extension =====
    rust_extensions=[
        RustExtension(
            "crisprme2._crisprme2_native",      # Python module path
            "native/rust/Cargo.toml",           # path to your Rust crate
            binding=Binding.PyO3,
            debug=False,
        ),
    ],
)
