# Releasing EnsembleQL

Releases are built in GitHub Actions. Maintainers do not upload distributions from a workstation and do not store PyPI API tokens in the repository.

## One-time account setup

1. Create `testpypi` and `pypi` environments under the GitHub repository's **Settings → Environments**. Require manual approval for `pypi`.
2. Register `ensembleql` as a pending Trusted Publisher on TestPyPI and PyPI with owner `eganeganegan`, repository `EnsembleQL`, workflow `release.yml`, and the corresponding environment name.
3. Connect the repository to Zenodo before the first GitHub release if a software DOI is desired.

## Prepare a release

1. Choose a version and update it in `pyproject.toml`, `CMakeLists.txt`, `python/ensembleql/__init__.py`, and `CITATION.cff`.
2. Move the version's changelog entry out of the unreleased section and add its release date.
3. Run `python tools/check_release.py` and the complete test suite.
4. Push the release-preparation commit and wait for every CI job to pass.

## Test the artifacts

Run the **Release** workflow manually with `publish_testpypi` disabled to build and test every artifact without publishing. When that succeeds, rerun it with `publish_testpypi` enabled. Trusted Publishing uploads the same source distribution and wheel matrix to TestPyPI.

TestPyPI and PyPI are separate indexes. A version uploaded to TestPyPI can still be uploaded to PyPI, but files on either index cannot be replaced. Increment the version if a published artifact needs to change.

## Publish

Create and push an annotated tag only after CI and the TestPyPI rehearsal pass:

```bash
git tag -a v0.1.0 -m "EnsembleQL 0.1.0"
git push origin v0.1.0
```

The tag must exactly match the version recorded in all public metadata. The workflow builds fresh artifacts, waits for approval in the `pypi` environment, publishes through Trusted Publishing, and creates a GitHub release with generated notes and attached distributions.

After publishing, install a wheel from PyPI in a clean environment, run a representative query, confirm the GitHub and Zenodo records, and start the next changelog section.
