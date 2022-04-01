#!/usr/bin/env python3
import logging
import os
import shutil
import re
import fire
from subprocess import run
from sh import sed

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("mvmapit-oscar-tasks")


def _cmd(command, args=[]):
    return run(command.split(" ") + args)


def _use_dev_files(tmpdir, files):
    logger.info(f"Make .dev files available: {files}")
    for f in files:
        src = os.path.join(tmpdir, f)
        dst = os.path.join(tmpdir, f.removesuffix(".dev"))
        shutil.copyfile(src, dst)
        logger.info(f"Copied {src} to {dst}")


def task_rinstall(tmpdir):
    logger.info(f"Temporary directory {tmpdir}.")
    dev_files = ["src/Makevars.in.dev", "src/MAPIT.h.dev"]
    _use_dev_files(tmpdir, dev_files)

    logger.info(f"Install mvMAPIT.")
    _cmd(f"R CMD INSTALL {tmpdir} --preclean")


def change_version(version):
    files = {
        "DESCRIPTION": f"s/(Version: )[0-9.]*/\\1{version}/",
        "configure.ac": f"s/[0-9]+\.[0-9]+\.[0-9]+/{version}/",
    }

    logger.info(f"Version: {version}")
    for key, sedcmd in files.items():
        sed(["-i", ".bu", "-E", sedcmd, key])

    _cmd("autoconf")


if __name__ == "__main__":
    fire.Fire()
