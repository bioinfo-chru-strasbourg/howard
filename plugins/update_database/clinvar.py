from howard.functions.commons import compress_file, command
from utils import (
    find_files,
    now,
    read_md5_file,
    recursive_chmod,
)
import tempfile
import logging as log
from os.path import join as osj
import os
import re
from plugins.update_database.factory import Database


class Clinvar(Database):
    def __init__(
        self,
        link=None,
        database=None,
        exclude_link=None,
        databases_folder=None,
        input=None,
        config_json=None,
        current_folder=None,
    ):
        super().__init__(
            link,
            database,
            exclude_link,
            databases_folder,
            input,
            config_json,
            current_folder,
        )

    def is_clinvar_up_to_date(
        self, assembly: str, clinvar_local_folder: str, clinvar_assembly: str
    ) -> bool:
        """
        Check if clinvar in local server is up to date, it compares the md5 in current version locally with the one on ucsc current version

        :param assembly: either hg19 or hg38
        :param clinvar_local_folder: current clinvar folder locally
        :param clinvar_assembly: either vcf_GRCh37 or vcf_GRCh38 on ucsc serveur
        :return bool: True or False regarding md5 differences
        """
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp_dir:
            log.info(f"Assembly {assembly}")
            clinvar_md5_link = osj(
                self.config_json["database"]["clinvar"],
                clinvar_assembly,
                "clinvar.vcf.gz.md5",
            )
            log.debug(clinvar_md5_link)
            clinvar_md5_local = osj(tmp_dir, "clinvar.vcf.gz.md5")
            clinvar_md5_current = find_files(
                osj(clinvar_local_folder, assembly), suffix=".md5"
            )[0]
            log.debug(f"Clinvar current: {clinvar_md5_current}")
            self.download(clinvar_md5_link, tmp_dir)
            if read_md5_file(clinvar_md5_local) != read_md5_file(clinvar_md5_current):
                log.warning(f"Clinvar not up to date {assembly}")
                return False
            else:
                log.info(
                    f"Clinvar version up to date: {' '.join(find_files(osj(clinvar_local_folder, assembly), suffix='md5'))}"
                )
                return True

    def get_clinvar_date(self):
        """
        Get date of current version of clinvar on ucsc server otherwise return None
        """
        list_files = self.list_databases(
            osj(self.config_json["database"]["clinvar"], "vcf_GRCh38")
        )
        for file in list_files:
            pattern = r"(?<=clinvar_)\d+(?=\.vcf\.gz)"
            match = re.search(pattern, file)
            if match:
                return match.group()
        return None

    def get_clinvar_assembly_folder(self, latest_folder_name: str, assembly: str):
        """
        Get assembly matching on ucsc server and associated path

        :param latest_folder_name: Name of new version of clinvar, should be a date
        :return clinvar_assembly: UCSC name of folder
        :return clinvar_assembly_folder: full path  locally of folder for a specific assembly
        """
        if assembly == "hg19":
            clinvar_assembly = "vcf_GRCh37"
        elif assembly == "hg38":
            clinvar_assembly = "vcf_GRCh38"
        clinvar_assembly_folder = osj(
            self.databases_folder, self.database, latest_folder_name, assembly
        )
        return clinvar_assembly, clinvar_assembly_folder

    def download_format_clinvar(
        self,
        clinvar_assembly,
        clinvar_assembly_folder,
        assembly,
        clinvar_local_folder,
    ):
        """
        Download and format clinvar vcf to parquet

        :param clinvar_assembly: UCSC name of folder
        :param clinvar_assembly_folder: full path  locally of folder for a specific assembly
        :param assembly: hg19 or hg38
        :param clinvar_local_folder: normally current folder for clinvar but could be specific
        """
        self.download_folder(
            self.list_databases(
                linkpage=osj(self.config_json["database"]["clinvar"], clinvar_assembly)
            ),
            clinvar_assembly_folder,
            osj(self.config_json["database"]["clinvar"], clinvar_assembly),
            starts="clinvar_",
            skip="papu",
        )
        # CHeckMD5
        if self.check_md5(
            find_files(osj(clinvar_local_folder, assembly), suffix=".md5")[0],
            find_files(osj(clinvar_local_folder, assembly), suffix=".vcf.gz")[0],
        ):
            self.formatting_clinvar(
                find_files(clinvar_assembly_folder, suffix=".vcf.gz")[0]
            )

    def update_clinvar_assembly(
        self, clinvar_local_folder, latest_folder_name, assembly
    ):
        """
        Update clinvar assembly only if md5 from ucsc latest version and current are different

        :param clinvar_local_folder: normally current folder for clinvar but could be specific
        :param latest_folder_name: Name of new version of clinvar, should be a date
        :param assembly: hg19 or hg38
        """
        # Both hg19 and hg38
        clinvar_assembly, clinvar_assembly_folder = self.get_clinvar_assembly_folder(
            latest_folder_name, assembly
        )
        if os.path.exists(clinvar_assembly_folder):
            log.info(
                f"Clinvar version {latest_folder_name} already exists locally EXIT"
            )
            exit()
        # Download folder
        if not self.is_clinvar_up_to_date(
            assembly, clinvar_local_folder, clinvar_assembly
        ):
            os.makedirs(clinvar_assembly_folder)
            self.download_format_clinvar(
                clinvar_assembly,
                clinvar_assembly_folder,
                assembly,
                clinvar_local_folder,
            )
        else:
            log.info(f"Clinvar version {latest_folder_name} is up to date EXIT")
            exit()
        # CHeckMD5
        if self.check_md5(
            find_files(osj(clinvar_local_folder, assembly), suffix=".md5")[0],
            find_files(osj(clinvar_local_folder, assembly), suffix=".vcf.gz")[0],
        ):
            self.formatting_clinvar(
                find_files(clinvar_assembly_folder, suffix=".vcf.gz")[0]
            )

    def update_clinvar(self):
        """
        Main process to update clinvar data
        """

        clinvar_local_folder = osj(
            self.databases_folder, self.database, self.current_folder
        )
        log.info(f"Clinvar folder {clinvar_local_folder}")
        if not clinvar_local_folder.endswith("current") and not os.path.exists(
            clinvar_local_folder
        ):
            raise (FileNotFoundError(f"{clinvar_local_folder} does not exist"))
        elif clinvar_local_folder.endswith("current") and not os.path.exists(
            clinvar_local_folder
        ):
            log.warning("Current clinvar folder does not exist, create latest version")

        latest_folder_name = self.get_clinvar_date()
        log.info(f"Clinvar current version on UCSC: {latest_folder_name}")
        if latest_folder_name is None:
            latest_folder_name = now()
            log.warning("Download clinvar, set version of the day")

        if not os.path.exists(clinvar_local_folder) and not os.path.exists(
            osj(self.databases_folder, self.database, latest_folder_name)
        ):
            for assembly in ["hg19", "hg38"]:
                clinvar_assembly, clinvar_assembly_folder = (
                    self.get_clinvar_assembly_folder(latest_folder_name, assembly)
                )
                self.download_format_clinvar(
                    clinvar_assembly,
                    clinvar_assembly_folder,
                    clinvar_assembly,
                    osj(
                        self.databases_folder,
                        self.database,
                        latest_folder_name,
                        assembly,
                    ),
                )
        else:
            for assembly in ["hg19", "hg38"]:
                self.update_clinvar_assembly(
                    clinvar_local_folder, latest_folder_name, assembly
                )

        # Update symlink
        log.info("Change permissions")
        recursive_chmod(
            osj(self.databases_folder, self.database, latest_folder_name), 0o755
        )
        log.info("Update symlink latest")
        try:
            os.unlink(osj(self.databases_folder, self.database, "latest"))
        except FileNotFoundError:
            log.debug("First latest symlink")
        os.symlink(
            latest_folder_name, osj(self.databases_folder, self.database, "latest")
        )
        log.warning(
            "You still need to validate the latest version of clinvar before using it in production"
        )

    def formatting_clinvar(self, clinvar_vcf_raw):
        """
        Add chr in front of contig name and transform resulting vcf into parquet, BGZIP required

        :param clinvar_vcf_raw: path of clinvar raw file from ucsc
        """
        clinvar_vcf = osj(os.path.dirname(clinvar_vcf_raw), "clinvar.vcf.gz")
        log.debug(f"Add 'chr' to {clinvar_vcf}")
        command(
            "zcat "
            + clinvar_vcf_raw
            + ' | awk \'{if ($0 !~ /^#/) $0="chr"$0; print}\' OFS="\t" > '
            + clinvar_vcf.replace(".vcf.gz", ".vcf")
        )
        compress_file(clinvar_vcf.replace(".vcf.gz", ".vcf"), clinvar_vcf)
        os.remove(clinvar_vcf.replace(".vcf.gz", ".vcf"))
        log.info(f"Formatting {clinvar_vcf} to parquet")
        self.vcf_to_parquet(clinvar_vcf)
