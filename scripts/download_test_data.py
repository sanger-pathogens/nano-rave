#!/usr/bin/env python3
# Copyright (C) 2022,2023 Genome Research Ltd.

import xml.etree.ElementTree as ET
import urllib.parse as urlparse
import urllib.request as urlrequest
from pathlib import Path
import shutil
from typing import Tuple, Optional, List, Any
import os


BUCKET_URL = "https://test-data-for-nano-rave.cog.sanger.ac.uk"
OUTPUT_DIR = Path("./test_data")


def _error_if_not_ok(response):
    status_code = response.getcode()
    # Note: 200 means OK
    if status_code != 200: 
        raise RuntimeError(f"Server responded with status code {status_code}")


def fetch_index(continuation_token = None) -> Tuple[List[str], Optional[Any]]:
    """
    Get the file path of every object in the bucket.

    Returns a tuple, the first element of which is the list of paths. The second
    is either None or a continuation token which means there are more objects still
    to be fetched. This continuation token can be passed back into this function to make
    a follow-up request to continue fetching objects.
    """

    url = BUCKET_URL + "?list-type=2"

    if continuation_token:
        url += "&continuation-token=" + continuation_token

    with urlrequest.urlopen(url) as response:
        _error_if_not_ok(response)
        etree = ET.parse(response)
    root = etree.getroot()

    # Get the particular xml namespace in use so we don't have to hardcode the version.
    # E.g. the namespace might be http://s3.amazonaws.com/doc/2006-03-01/, in which case
    # the root tag will be "{http://s3.amazonaws.com/doc/2006-03-01}ListBucketResult"
    s3_namespace = root.tag.split("}")[0].strip("{")
    ns = { "s3": s3_namespace }

    contents = root.findall("s3:Contents", ns)
    file_paths = [object.find("s3:Key", ns).text for object in contents]
    
    is_truncated = root.find("s3:IsTruncated", ns).text == "true"
    if is_truncated:
        continuation_token = root.find("s3:NextContinuationToken", ns).text
        return (file_paths, continuation_token)
    else:
        return (file_paths, None)


def download_file(path: str):
    url = BUCKET_URL + "/" + urlparse.quote_plus(path)

    # Note: '/' operator joins paths
    dest = OUTPUT_DIR / path
    os.makedirs(dest.parent, exist_ok=True)
    
    with urlrequest.urlopen(url) as response, open(dest, 'wb') as out_file:
        _error_if_not_ok(response)
        shutil.copyfileobj(response, out_file)


if __name__ == "__main__":
    if os.path.exists(OUTPUT_DIR):
        print(f"Error: output directory {OUTPUT_DIR} already exists")
        exit(1)

    continuation_token = None
    while True:
        print("Fetching ", end="")
        if continuation_token:
            print("continued ", end="")
        print(f"object list from {BUCKET_URL}...", end="")
        
        (paths, continuation_token) = fetch_index(continuation_token)
        print(f" Got {len(paths)} objects")
        for path in paths:
            print(f"Downloading {path} into {OUTPUT_DIR}...", end="")
            download_file(path)
            print(" Done")
        
        if continuation_token is None:
            break
    
    print(f"All test data downloaded into {OUTPUT_DIR}")
