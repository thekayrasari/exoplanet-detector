# Contributing to Exoplanet Detector

Thank you for considering contributing! Here is how to get involved.

## 🐛 Reporting Bugs

Open an issue and include:
- Your operating system and Python version
- The full error message / traceback
- The target name that caused the problem
- Steps to reproduce

## 💡 Suggesting Features

Open an issue with the label **enhancement** and describe:
- What you want to do that the tool cannot do today
- Why it would be useful to other users

## 🔧 Submitting Code

1. Fork the repository and create a branch from `master`
2. Keep all code inside `exoplanet_pipeline.py` (single-file design)
3. Run the self-test before opening a PR:
   ```bash
   python exoplanet_pipeline.py self-test
   ```
4. Write a clear commit message describing what changed and why
5. Open a Pull Request — small, focused changes are easier to review

## 📐 Code Style

- PEP 8 formatting
- Docstrings on all public methods
- Prefer clarity over cleverness — this is a science tool used by non-programmers

## 🔭 Science Contributions

If you discover a new exoplanet candidate using this tool, please:
1. Report it to [NASA ExoFOP](https://exofop.ipac.caltech.edu/tess/) for follow-up
2. Feel free to open an issue here linking to your ExoFOP submission!

## 📄 License

By contributing, you agree that your contributions will be licensed under the MIT License.
