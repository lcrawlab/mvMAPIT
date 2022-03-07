#!/usr/bin/env python3
import logging
import os
import shutil
from subprocess import run

import git
from jinja2 import Template

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("mvmapit-oscar-tasks")


def _cmd(command, args=[]):
    return run(command.split(" ") + args)


def task_rinstall(tmpdir):
    logger.info(f"Temporary directory {tmpdir}.")
    logger.info(f"Install mvMAPIT.")
    _cmd(f"R CMD INSTALL {tmpdir} --preclean")
