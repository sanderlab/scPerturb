
from pathlib import Path
from typing import Union
from urllib.parse import (urlparse, unquote)
import subprocess
import requests
import shlex


# addresses #
#############

URLS = {
    "chrom_modifiers": [
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4887nnn/GSM4887677/suppl/GSM4887677_screen1_snATAC.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4887nnn/GSM4887678/suppl/GSM4887678_screen1_snATACguide.IDs.mat.txt.gz"
    ],
    "chrom_remodelers": [
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4887nnn/GSM4887679/suppl/GSM4887679_screen2_snATAC.bed.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4887nnn/GSM4887680/suppl/GSM4887680_screen2_snATACguide.IDs.mat.txt.gz"
    ]
}




def download_file(
    url: str,
    output_dir: Union[Path, str],
    overwrite: bool = False,
    create_new_file: bool = False,
) -> Path:
    """Downloads the content behind the specified url and writes the bytes of the content into a file
    inside the specified output directory.
    The filename is determined from the url instead of the url header's content-disposition, because it is
    not always available. This also ensures that the function gives the same filename as a call to wget:

    The Content-Disposition header can be used by a server to suggest a filename for a downloaded file.
    By default, wget uses the last part of the URL as the filename,
    but you can override this with --content-disposition,
    which uses the server's suggested name.

    :param url: url to downloadable content
    :type url: str
    :param output_dir: Path to output directory.
    :type output_dir: Union[Path,str]
    :param overwrite: Determines whether an existing file will be overwritten.
    If False, will add a numerical suffix to the filename.
    :type overwrite: bool
    :param create_new_file: Only relevant if overwrite is False.
    If True it will add a numerical suffix to the filename. Otherwise it will not write the file.
    :type create_new_file: bool
    :return: The file path to the downloaded file.
    :rtype: Path
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # download content
    req = requests.get(url)

    filename = urlparse(req.url).path.split("/")[-1]
    filename = unquote(filename)
    filepath = output_dir / filename

    if not overwrite:
        if create_new_file:
            # check if file exists and if so find a new filename.
            suffix_int = 1
            while filepath.is_file():
                filename_split = filename.rsplit(maxsplit=1)
                filename = f"{filename_split[0]}_{suffix_int}{filename_split[1]}"
                filepath = output_dir / filename
                suffix_int += 1
        else:
            if filepath.is_file():
                return filepath
    # write content to file
    filepath.write_bytes(req.content)
    return filepath

def run():
    output_dir = Path(__file__).parent / "data" / "CRISPR_sciATAC"
    for urls in URLS.values():
        for url in urls:
            filepath = download_file(url, output_dir, overwrite=False, create_new_file=False)

            if url.endswith('.bed.gz'):
                #the CRISPR_sciATAC files are not in tabix format. So annoying. I need them to be in tabix format for ArchR
                subprocess.call(['gunzip', str(filepath)])
                filepath = filepath.with_suffix('') #remove .gz
                sort_command = shlex.split(f'sort -k1,1V -k2,2n -k3,3n {filepath}')
                ps = subprocess.run(sort_command, check=True, capture_output=True)
                with open(filepath, "wb") as sorted_fragments:
                    sorted_fragments.write(ps.stdout)
                subprocess.call(['bgzip',  str(filepath)]) #special compression from samtools
                filepath = filepath.parent / (filepath.name + '.gz') #bgzip adds the gz suffix
                subprocess.call(['tabix',  str(filepath)]) #finally create index file with tabix


if __name__ == '__main__':
    run()