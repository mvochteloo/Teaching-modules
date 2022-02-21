#!/usr/bin/env python3

"""
File:         parse_results.py
Created:      2021/02/21
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
from __future__ import print_function
import argparse
import os

# Third party imports.
import pandas as pd

# Local application imports.

# Metadata
__program__ = "Parse Results"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.data_path = getattr(arguments, 'data')

        if not os.path.exists(self.data_path):
            print("Input file does not exist.")
            exit()
        if not self.data_path.endswith(".csv"):
            print("Input file should have the '.csv' extension.")
            exit()

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit.")
        parser.add_argument("-d",
                            "--data",
                            type=str,
                            required=True,
                            help="The path to the data file (.xlsx).")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        df = self.load_file(self.data_path, header=0, index_col=None)

        print("")
        print("### Step2 ###")
        print("Pre-processing data")
        preference_df = self.preprocess(df=df)

        print("")
        print("### Step3 ###")
        print("Saving data")
        self.save_file(df=preference_df,
                       outpath=self.data_path.replace(".csv", ".xlsx"))

    @staticmethod
    def load_file(path, header, index_col, sep=","):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def preprocess(df):
        data = []
        index = []
        columns = []
        for _, row in df.iterrows():
            student_id = "{}, {} [{}]".format(row["Last Name"], row["First Name"], str(row["Username"]))
            preference = {}
            for score, answer in enumerate(row["Answer"].split(",")):
                project_id = answer.split("]")[0].replace("[", "")
                preference[project_id] = (score + 1)

            student_answer = []
            for i in range(1, 15):
                project_id = "PROJECT{}".format(i)
                if project_id in preference.keys():
                    if project_id not in columns:
                        columns.append(project_id)
                    student_answer.append(preference[project_id])

            data.append(student_answer)
            index.append(student_id)

        preference_df = pd.DataFrame(data, index=index, columns=columns)
        preference_df.sort_index(inplace=True)

        return preference_df

    @staticmethod
    def save_file(df, outpath, header=True, index=True, na_rep="",
                  sheet_name="Sheet1"):
        df.to_excel(outpath,
                    sheet_name=sheet_name,
                    na_rep=na_rep,
                    header=header,
                    index=index)

        print("\tSaved dataframe: {} "
              "with shape: {}".format(os.path.basename(outpath),
                                      df.shape))

    def print_arguments(self):
        print("Arguments:")
        print("  > Data file: {}".format(self.data_path))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
