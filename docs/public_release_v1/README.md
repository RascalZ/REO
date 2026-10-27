# Version 0.1.0 publication evidence

The current machine-readable publication state is
[PUBLICATION_STATUS.json](PUBLICATION_STATUS.json). The GitHub repository is
public. The [GitHub release](https://github.com/RascalZ/REO/releases/tag/v0.1.0) and
[Zenodo archive](https://doi.org/10.5281/zenodo.23681422) are published.
Anonymous downloads are verified separately from publication metadata.

- [Payload scope and checksums](PAYLOAD_SUMMARY.json)
- [Per-file catalog](PAYLOAD_ASSET_CATALOG.json)
- [Rebuild the model and metadata packages](PAYLOAD_BUILD.md)
- [Release notes](RELEASE_NOTES.md)
- [Asset license scope](MODEL_ASSETS_LICENSE.md)
- [Third-party attribution](THIRD_PARTY_ATTRIBUTION.md)

The release preserves 64 trained model payloads and four original fold files.
Four minimal identity/membership projections retain all rows and selected
fields but omit source experiment, curve, QC and dose pointers. This permits
fixed-model inference and assignment auditing; it does not supply raw training
data or the complete historical external-rescoring inputs. Upstream data
rights and compatibility dependencies remain separate.

The v4 author manuscripts and submission ZIP remain unchanged. The archival
DOI identifies software and research assets, not the unpublished paper. Once
the real release is verified, its access wording can replace the older
private/pending wording in a separately reviewed manuscript copy.

All publication work is limited to packaging, provenance, licensing,
synthetic-input software checks and download verification. No scientific
training, tuning, resampling or external drug-response scoring is repeated.
