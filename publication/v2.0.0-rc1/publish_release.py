#!/usr/bin/env python3
"""Publish only a fully verified, version-bound CosmicScanner pre-release."""
import argparse
import hashlib
import json
import os
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener


ROOT = Path(__file__).resolve().parent
REPOSITORY = "huangjie2025/CosmicScanner"
API = "https://api.github.com/repos/" + REPOSITORY


class SafeRedirects(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected and urlparse(req.full_url).netloc != urlparse(newurl).netloc:
            redirected.remove_header("Authorization")
        return redirected


def local_checks():
    manifest = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["tag"] == "v2.0.0-rc1"
    assert manifest["source_commit"] == "e80062ae31942dc518b8a13110c509470bfe7421"
    assert manifest["prerelease"] is True and manifest["make_latest"] is False
    names = [asset["name"] for asset in manifest["assets"]]
    assert len(names) == len(set(names)) == 10
    assert all(Path(name).name == name for name in names)
    for asset in manifest["assets"]:
        path = ROOT / "assets" / asset["name"]
        assert path.stat().st_size == asset["bytes"], path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == asset["sha256"], path
    notes = (ROOT / "RELEASE_NOTES.md").read_text(encoding="utf-8")
    assert "pre-release" in notes and "experimental" in notes
    return manifest, notes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    manifest, notes = local_checks()
    print("Local release verification: PASS (10 exact assets)", flush=True)
    if args.check_only:
        return
    assert os.environ.get("GITHUB_REPOSITORY") == REPOSITORY
    assert os.environ.get("GITHUB_REF") == "refs/heads/codex/v2.0.0-rc1-publication-20261020"
    token = os.environ["GITHUB_TOKEN"]
    opener = build_opener(SafeRedirects())

    def request(url, method="GET", payload=None, binary=None, accept=None, content_type=None):
        headers = {"Accept": accept or "application/vnd.github+json",
                   "Authorization": "Bearer " + token,
                   "X-GitHub-Api-Version": "2026-03-10",
                   "User-Agent": "CosmicScanner-release-verification"}
        data = binary
        if payload is not None:
            data = json.dumps(payload).encode("utf-8")
            headers["Content-Type"] = "application/json"
        elif binary is not None:
            headers["Content-Type"] = content_type or "application/octet-stream"
        req = Request(url, data=data, headers=headers, method=method)
        with opener.open(req, timeout=180) as response:
            raw = response.read()
        return raw if accept == "application/octet-stream" else json.loads(raw)

    tag = manifest["tag"]
    ref_url = API + "/git/ref/tags/" + quote(tag, safe="")
    try:
        tag_object = request(ref_url)
    except HTTPError as error:
        if error.code != 404:
            raise
        request(API + "/git/refs", "POST", {"ref": "refs/tags/" + tag,
                                            "sha": manifest["source_commit"]})
        tag_object = request(ref_url)
    assert tag_object["object"]["type"] == "commit", "Unexpected existing annotated tag; inspect before publishing"
    assert tag_object["object"]["sha"] == manifest["source_commit"], "Existing tag has a different source revision"

    releases = request(API + "/releases?per_page=100")
    release = next((entry for entry in releases if entry["tag_name"] == tag), None)
    if release is None:
        release = request(API + "/releases", "POST", {
            "tag_name": tag, "target_commitish": manifest["source_commit"],
            "name": manifest["title"], "body": notes, "draft": True,
            "prerelease": True, "make_latest": "false", "generate_release_notes": False})
    assert release["prerelease"] is True, "Existing release is not a pre-release"
    assert not release.get("immutable"), "Existing release is immutable"
    upload_url = release["upload_url"].split("{", 1)[0]
    assert urlparse(upload_url).netloc == "uploads.github.com"
    existing = request(release["assets_url"] + "?per_page=100")
    expected_names = {entry["name"] for entry in manifest["assets"]}
    assert {entry["name"] for entry in existing} <= expected_names, "Unexpected existing release assets"
    by_name = {entry["name"]: entry for entry in existing}
    verified = []
    for expected in manifest["assets"]:
        name = expected["name"]
        asset = by_name.get(name)
        if asset is None:
            print("Uploading " + name, flush=True)
            asset = request(upload_url + "?name=" + quote(name, safe=""), "POST",
                            binary=(ROOT / "assets" / name).read_bytes())
        assert asset["size"] == expected["bytes"], name
        assert asset["state"] == "uploaded", name
        digest = asset.get("digest")
        if digest:
            assert digest == "sha256:" + expected["sha256"], name
        downloaded = request(asset["url"], accept="application/octet-stream")
        assert len(downloaded) == expected["bytes"], name
        assert hashlib.sha256(downloaded).hexdigest() == expected["sha256"], name
        verified.append({"name": name, "sha256": expected["sha256"],
                         "bytes": len(downloaded), "id": asset["id"],
                         "url": asset["browser_download_url"]})
        print("Downloaded SHA-256 verified: " + name, flush=True)

    assert request(ref_url)["object"]["sha"] == manifest["source_commit"]
    published = request(release["url"], "PATCH", {
        "name": manifest["title"], "body": notes, "draft": False,
        "prerelease": True, "make_latest": "false"})
    assert published["draft"] is False and published["prerelease"] is True
    assert {entry["name"] for entry in published["assets"]} == expected_names
    receipt = {"status": "PASS", "repository": REPOSITORY, "tag": tag,
               "source_commit": manifest["source_commit"],
               "release_id": published["id"], "url": published["html_url"],
               "draft": published["draft"], "prerelease": published["prerelease"],
               "verification": "All 10 uploaded assets downloaded and SHA-256 matched before publication",
               "assets": verified}
    (ROOT / "published-release.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
