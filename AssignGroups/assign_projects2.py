#!/usr/bin/env python3

"""
File:         assign_projects.py
Created:      2021/02/02
Last Changed: 2022/02/21
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
import itertools
import os
import math

# Third party imports.
import pandas as pd
import seaborn as sns
from scipy import optimize

# Local application imports.

# Metadata
__program__ = "Assign Projects"
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
        self.sheet_name = getattr(arguments, 'sheet')
        self.min_group_size = getattr(arguments, 'min_group_size')
        self.max_group_size = getattr(arguments, 'max_group_size')

        if not os.path.exists(self.data_path):
            print("Input file does not exist.")
            exit()
        if not self.data_path.endswith(".xlsx"):
            print("Input file should be an excel sheet.")
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
        parser.add_argument("-s",
                            "--sheet",
                            type=str,
                            required=False,
                            default="Sheet1",
                            help="The name of the excel sheet. "
                                 "Default: 'Sheet1'.")
        parser.add_argument("-min",
                            "--min_group_size",
                            type=int,
                            required=False,
                            default=3,
                            help="The minimal number of student in a"
                                 "group. Default: 3.")
        parser.add_argument("-max",
                            "--max_group_size",
                            type=int,
                            required=False,
                            default=4,
                            help="The maximal number of student in a"
                                 "group. Default: 4.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        preference_df = self.load_file(self.data_path,
                                       header=0,
                                       index_col=0,
                                       sheet_name=self.sheet_name)

        print("")
        print("### Step2 ###")
        print("Validate data")
        self.validate(preference_df)

        print("")
        print("### Step3 ###")
        print("Pre-process the data")
        preference_df = preference_df.astype(int)
        print(preference_df)

        print("")
        print("### Step4 ###")
        print("Calculate information")
        students = list(preference_df.index)
        print("\tNumber of students: {}".format(len(students)))
        projects = list(preference_df.columns)
        print("\tNumber of projects: {}".format(len(projects)))
        identical = self.find_identical_answers(preference_df)
        print("\tStudents with identical answers:")
        if len(identical) > 0:
            for identical_answer_group in identical:
                print("\t\t{}".format(', '.join(identical_answer_group)))
        else:
            print("\t\tNone")
        avg_preference_df = preference_df.mean(axis=0)
        avg_preference_df.sort_values(inplace=True)
        print("Average score:")
        for project_id, mean_score in avg_preference_df.iteritems():
            print("\t{}: {:.1f}".format(project_id, mean_score))
        self.plot_preferences(preference_df)

        print("")
        print("### Step5 ###")
        print("Calulate group sizes")
        possible_groups_sizes = self.create_group_sizes(students, self.min_group_size, self.max_group_size)
        for pgs in possible_groups_sizes:
            counts = pd.Series(pgs).value_counts()
            text = "\t{} groups:\t".format(len(pgs))
            for index, value in counts.iteritems():
                text += "{} with {} students\t".format(value, index)
            print(text)

        print("")
        print("### Step6 ###")
        print("Solve fast")
        solutions = self.solve(preference_df, possible_groups_sizes)
        for i, (sum_score, groups) in enumerate(solutions):
            print("Solution {}:".format(i + 1))

            total_choice_counts = {}
            for project, students in groups.items():
                choice_counts = {}
                students_str = []
                for name, score in students:
                    students_str.append("{} [{}]".format(name, score))
                    if score in choice_counts:
                        choice_counts[score] = choice_counts[score] + 1
                    else:
                        choice_counts[score] = 1

                if len(students_str) > 0:
                    choice_info = []
                    for j in range(1, 100):
                        if j in choice_counts:
                            choice_info.append("{}x{}".format(j, choice_counts[j]))
                    print("\t  {}: {}, ({})".format(project, ", ".join(students_str), ", ".join(choice_info)))
                else:
                    print("\t  {}".format(project))

                for choice, count in choice_counts.items():
                    if choice in total_choice_counts:
                        total_choice_counts[choice] += count
                    else:
                        total_choice_counts[choice] = count

            choice_info = []
            for j in range(1, 100):
                if j in total_choice_counts:
                    choice_info.append("{}x{}".format(j, total_choice_counts[j]))

            print("\t  Preference score: {}\t\tChoice counts: {}".format(sum_score, ", ".join(choice_info)))

    @staticmethod
    def load_file(path, header=None, index_col=None, nrows=None, skiprows=0, sheet_name=None):
        df = pd.read_excel(path, header=header, index_col=index_col,
                           nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        return df

    @staticmethod
    def validate(df):
        if len(df.index.unique()) != df.shape[0]:
            print("\tThe 'Student' column cannot duplicate values.")
            exit()

        for i, (student_id, student_preferences) in enumerate(df.iterrows()):
            preferences_counts = student_preferences.value_counts()
            preferences_counts.sort_index(inplace=True)
            prev_rank = 1
            for rank, count in preferences_counts.iteritems():
                if rank > (prev_rank + 1):
                    print("\t[{}] {} does not have a valid top preference".format(i, student_id))
                    exit()
                prev_rank = rank

            # student_preferences = list(zip(student_preferences.index, student_preferences.values))
            # student_preferences.sort(key=lambda x: x[1])
            #
            # if student_id == ":
            #     print(student_id)
            #     for project, score in student_preferences:
            #         print("{}. [{}] ...".format(score, project))
            #     print("")

    @staticmethod
    def find_identical_answers(df):
        choice = {}
        for student, row in df.iterrows():
            order = "".join([str(x) for x in row.values])
            if order in choice.keys():
                students_with_this_choice = choice[order]
                students_with_this_choice.append(student)
                choice[order] = students_with_this_choice
            else:
                choice[order] = [student]

        # filter.
        students_with_identical_answers = []
        for item in choice.values():
            if len(item) > 1:
                students_with_identical_answers.append(item)
        return students_with_identical_answers

    @staticmethod
    def plot_preferences(df):
        counts = df.apply(pd.value_counts, axis=0)
        counts.columns = [x[:20] for x in counts.columns]
        counts.fillna(0, inplace=True)
        counts.reset_index(drop=False, inplace=True)
        counts = counts.melt(id_vars=["index"])

        sns.set_theme(style="whitegrid")
        g = sns.catplot(
            data=counts, kind="bar", orient="h", ci=None,
            x="value", y="index", hue="variable", aspect=.6
        )
        g.despine(left=True)
        g.set_axis_labels("count", "preference")
        g.legend.set_title("")
        g.savefig("preferences.png")

    @staticmethod
    def create_group_sizes(students, min_group_size, max_group_size):
        possible_groups_sizes = []
        for group_size in range(min_group_size, max_group_size + 1, 1):
            n_groups = math.floor(len(students) / group_size)
            groups = [group_size for _ in range(n_groups)]
            students_left = len(students) - (n_groups * group_size)

            # option 1: divide students of the groups.
            if students_left < len(groups):
                tmp_groups = groups.copy()
                for i, student in enumerate(range(students_left)):
                    tmp_groups[i] += 1
                if max(tmp_groups) <= max_group_size and tmp_groups not in possible_groups_sizes:
                    possible_groups_sizes.append(tmp_groups)

            # option 2: add an extra group with the remainder.
            if students_left >= min_group_size:
                tmp_groups = groups.copy()
                tmp_groups.append(students_left)
                if min(tmp_groups) >= min_group_size and tmp_groups not in possible_groups_sizes:
                    possible_groups_sizes.append(tmp_groups)

        return possible_groups_sizes

    @staticmethod
    def solve(preference, valid_group_sizes):
        n_projects = preference.shape[1]

        solutions = []
        for group_sizes in valid_group_sizes:
            # Make sure the group sizes equals the number of projects.
            sizes = group_sizes.copy()
            n_filled_projects = len(sizes)
            for _ in range(n_projects - n_filled_projects):
                sizes.append(None)

            # Get all possible combinations of project - group size.
            permutations = list(itertools.permutations(sizes))
            found = set()
            project_sizes = []
            for permutation in permutations:
                permutation_str = "".join([str(x) if x is not None else " " for x in permutation])
                if permutation_str in found:
                    continue

                project_sizes.append(permutation)
                found.add(permutation_str)

            # Loop through the options.
            for project_size in project_sizes:

                # Construct the cost matrix.
                cost_list = []
                for column_index, repeats in enumerate(project_size):
                    if repeats is None:
                        continue
                    cost_df = preference.iloc[:, [column_index] * repeats]
                    cost_df.columns = ["{} - pos {}".format(col, i) for i, col in enumerate(cost_df.columns)]
                    cost_list.append(cost_df)
                cost_df = pd.concat(cost_list, axis=1)

                # Solve.
                solution = optimize.linear_sum_assignment(cost_df)

                # Construct the groups.
                groups = {project: [] for project in preference.columns}
                sum_score = 0
                for row_id, col_id in zip(*solution):
                    student_id = cost_df.index[row_id]
                    project_id = cost_df.columns[col_id].split(" -")[0]
                    score = preference.loc[student_id, project_id]
                    groups[project_id].append((student_id, score))
                    sum_score += score

                # Save.
                solutions.append([sum_score, groups])

        # Sort the solutions.
        solutions.sort(key=lambda x: x[0])

        return solutions

    def print_arguments(self):
        print("Arguments:")
        print("  > Data file: {}".format(self.data_path))
        print("  > Sheet name: {}".format(self.sheet_name))
        print("  > Min group size: {}".format(self.min_group_size))
        print("  > Max group size: {}".format(self.max_group_size))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
