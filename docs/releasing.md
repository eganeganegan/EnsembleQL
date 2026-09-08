# Releasing EnsembleQL

Releases are built in GitHub Actions. Maintainers do not upload distributions from a workstation and do not store PyPI API tokens in the repository.

## One-time account setup

1. Create `testpypi` and `pypi` environments under the GitHub repository's **Settings → Environments**. Require manual approval for `pypi`.
2. Register `ensembleql` as a pending Trusted Publisher on TestPyPI and PyPI with owner `eganeganegan`, repository `EnsembleQL`, workflow `release.yml`, and the corresponding environment name.
3. Connect the repository to Zenodo before the first GitHub release if a software DOI is desired. Sign in to Zenodo with the maintainer account, open the profile menu's **GitHub** page, click **Sync now**, find `eganeganegan/EnsembleQL`, and enable its toggle. Open the repository in Zenodo and select **Create release** if that action is shown.

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

## Zenodo archive and DOI

The production tag workflow creates the GitHub release after PyPI accepts every artifact. When the repository has been enabled in Zenodo beforehand, Zenodo ingests that GitHub release and creates a versioned software record. Processing can take a few minutes.

After the record appears:

1. Check its title, authors, license, version, description, and release date. Zenodo reads the repository's validated `CITATION.cff` metadata.
2. Confirm that the uploaded archive corresponds to the `v0.1.0` GitHub release.
3. Copy the DOI into papers and software metadata that require a specific release. Prefer the concept DOI when citing EnsembleQL generally and the version DOI when exact reproducibility matters.
4. On the Zenodo record, select **Get the DOI badge**, copy its Markdown, and replace the temporary `DOI: pending` badge at the top of `README.md`. Use the concept DOI badge if it should continue representing all EnsembleQL releases.
