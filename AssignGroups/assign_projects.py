#!/usr/bin/env python3

"""
File:         assign_projects.py
Created:      2021/02/02
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
import itertools
import os
import math
import multiprocessing as mp
import queue
import random
import time
from datetime import datetime

# Third party imports.
import numpy as np
import pandas as pd
import cvxpy as cp
import seaborn as sns

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
        cores = getattr(arguments, 'cores')
        if cores is None:
            cores = mp.cpu_count()
        self.cores = cores
        self.max_end_time = int(time.time()) + self.time_to_sec(getattr(arguments, 'time'))

        if not os.path.exists(self.data_path):
            print("Input file does not exist.")
            exit()
        if not self.data_path.endswith(".xlsx"):
            print("Input file should be an excel sheet.")
            exit()

    def time_to_sec(self, time):
        if time == "SHORT":
            return 21300
        elif time == "MEDIUM":
            return 86100
        elif time == "LONG":
            return 604500
        else:
            print("Unexpected input for -t / --time.")
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
        parser.add_argument("-c",
                            "--cores",
                            type=int,
                            required=False,
                            default=None,
                            help="The number of cores to use. Default: max.")
        parser.add_argument("-t",
                            "--time",
                            type=str,
                            required=False,
                            choices=["SHORT", "MEDIUM", "LONG"],
                            default="SHORT",
                            help="The maximal time for the program to run. "
                                 "Default: 6 hours.")

        return parser.parse_args()

    def start(self):
        self.print_arguments()

        print("### Step1 ###")
        print("Loading data")
        df = self.load_file(self.data_path, header=0, sheet_name=self.sheet_name)

        print("")
        print("### Step2 ###")
        print("Validate data")
        self.validate(df)

        print("")
        print("### Step3 ###")
        print("Split data")
        df.set_index("Student", inplace=True)
        info = df[["Bachelor"]].copy()
        print(info)
        preference = df.iloc[:, 1:].copy()
        print(preference)
        del df

        print("")
        print("### Step4 ###")
        print("Calculate information")
        students = list(info.index)
        print("\tNumber of students: {}".format(len(students)))
        bachelors = list(info.loc[:, "Bachelor"].unique())
        print("\tNumber of unique bachelor studies: {}".format(len(bachelors)))
        projects = list(preference.columns)
        print("\tNumber of projects: {}".format(len(projects)))
        student_to_study_dict = dict(zip(info.index, info.loc[:, "Bachelor"]))
        identical = self.find_identical_answers(preference)
        print("\tStudents with identical answers:")
        if len(identical) > 0:
            for identical_answer_group in identical:
                print("\t\t{}".format(', '.join(identical_answer_group)))
        else:
            print("\t\tNone")
        # self.print_preferences(preference)
        self.plot_preferences(preference)

        print("")
        print("### Step5 ###")
        print("Calulate group sizes")
        possible_groups_sizes = self.create_group_sizes(students, self.min_group_size, self.max_group_size)
        partitions = []
        for pgs in possible_groups_sizes:
            partitions.append(len(pgs))
            counts = pd.Series(pgs).value_counts()
            text = "\t{} groups:\t".format(len(pgs))
            for index, value in counts.iteritems():
                text += "{} with {} students\t".format(value, index)
            print(text)

        print("")
        print("### Step6 ###")
        print("Solve fast")
        solutions = self.solve(preference, possible_groups_sizes)
        for i, solution in enumerate(solutions):
            print("Solution {}:".format(i))
            preference_score = 0
            diversity_score = 0
            for project, members in solution.items():
                text = []
                studies = {}
                dev_score = 0
                for member in members:
                    study = student_to_study_dict[member]
                    if study in studies:
                        n = studies[study]
                        studies[study] = n + 1
                    else:
                        dev_score += 1
                        studies[study] = 1
                    pref_score = preference.at[member, project]
                    preference_score += pref_score
                    text.append("{} [{} choice]".format(member, pref_score))
                diversity_score += dev_score
                print("\t  {}: {} {}".format(project, ", ".join(text), studies))
            print("\t  Preference score: {}\t\tDiversity score: {}".format(preference_score, diversity_score))

        print("")
        print("### Step7 ###")
        print("Create unique groups")
        possible_groups = []
        options = 0
        for k in partitions:
            for groups in self.sorted_k_partitions(students, k):
                group_size = [len(x) for x in groups]
                if group_size in possible_groups_sizes:
                    possible_groups.append(groups)
                    options += len([list(zip(x, groups)) for x in itertools.permutations(preference.columns, len(groups))])
                    # sorted_groups = [sorted(x) for x in groups]
                    # sorted_groups.sort(key=lambda x: x[1])
                    # if sorted_groups not in possible_groups:
                    #     possible_groups.append(sorted_groups)
        print("\tThere are {} possible groups arrangements.".format(len(possible_groups)))
        print("\tThere are {} possible group - project arrangements.".format(options))

        print("")
        print("### Step8 ###")
        print("Score groups")
        my_mp = MultiProcessor(cores=self.cores,
                               max_end_time=self.max_end_time,
                               data=[info, preference])
        scores = my_mp.process(self.score_group, possible_groups)

        print("Sorting {} groups.".format(len(scores)))
        scores.sort(key=lambda x: (x[1], -x[2]))

        print("")
        print("### Step8 ###")
        print("Finding the best score per criteria.")
        best_pref_score = np.inf
        prev_pref_score = None
        best_div_score = -np.inf
        for i, (group_assignment, preference_score, diversity_score) in enumerate(scores):
            if (preference_score <= best_pref_score) or (preference_score <= best_pref_score and diversity_score > best_div_score):
                print("Solution {}:".format(i))
                group_assignment.sort(key=lambda x: x[0])

                choice_counts = {}

                if preference_score < best_pref_score:
                    best_pref_score = preference_score
                if diversity_score > best_div_score:
                    best_div_score = diversity_score

                for project, members in group_assignment:
                    text = []
                    studies = {}
                    dev_score = 0
                    for member in members:
                        study = student_to_study_dict[member]
                        if study in studies:
                            n = studies[study]
                            studies[study] = n + 1
                        else:
                            dev_score += 1
                            studies[study] = 1
                        pref_score = preference.at[member, project]
                        if pref_score in choice_counts:
                            choice_counts[pref_score] = choice_counts[pref_score] + 1
                        else:
                            choice_counts[pref_score] = 1

                        text.append("{} [{} choice]".format(member, pref_score))
                    print("\t  {}: {} {}".format(project, ", ".join(text), studies))
                choice_info = []
                for j in range(1, 100):
                    if j in choice_counts:
                        choice_info.append("{}x{}".format(j, choice_counts[j]))

                print("\t  Preference score: {}\t\tDiversity score: {}\t\tN Groups: {}\t\tChoice counts: {}".format(preference_score, diversity_score, len(group_assignment), ", ".join(choice_info)))
                i += 1

            if prev_pref_score is not None and prev_pref_score < preference_score:
                break

            prev_pref_score = preference_score

        print("")
        print("### Step9 ###")
        print("Saving the top 50 options.")
        good_options = pd.DataFrame(scores[:50], columns=["group_assignment", "preference_score", "diversity_score"])
        good_options.to_excel("good_group_distributions.xlsx")

    @staticmethod
    def load_file(path, header=None, index_col=None, nrows=None, skiprows=0, sheet_name=None):
        df = pd.read_excel(path, header=header, index_col=index_col,
                           nrows=nrows, skiprows=skiprows, sheet_name=sheet_name)
        return df

    @staticmethod
    def validate(df):
        for (index, colname, description, canHaveDuplicates) in [(0, "Student", "student name", False),
                                                                 (1, "Bachelor", "bachelor's degree", True)]:
            if df.columns[index] != colname:
                print("\tThe first column should be named '{}' and contain "
                      "the {}.".format(colname, description))
                exit()
            if df.iloc[:, index].isnull().any():
                print("\tThe '{}' column cannot contain missing values.".format(colname))
                exit()
            if not canHaveDuplicates:
                if len(df.iloc[:, index].unique()) != df.shape[0]:
                    print("\tThe '{}' column cannot duplicate values.".format(colname))
                    exit()

        for i, (_, row) in enumerate(df.iterrows()):
            counts = row[2:].value_counts()
            counts.sort_index(inplace=True)
            prev_rank = 1
            for rank, count in counts.iteritems():
                if rank > (prev_rank + 1):
                    print("\t[{}] {} - {} does not have a valid top three".format(i, row["Student"], row["Bachelor"]))
                    exit()
                prev_rank = rank

    @staticmethod
    def find_identical_answers(preference):
        choice = {}
        for student, row in preference.iterrows():
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
    def print_preferences(df):
        tmp = df.copy()
        for index, row in tmp.iterrows():
            print(index)
            row.sort_values(inplace=True)
            print(row)
            print("")

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
    def construct_background_matrix(df):
        bg_df = pd.DataFrame(0, index=df.index, columns=df.loc[:, "Bachelor"].unique())
        for student in df.index:
            bg_df.at[student, df.at[student, "Bachelor"]] = 1
        return bg_df

    @staticmethod
    def solve(preference, valid_group_sizes):
        solutions = []

        selection = cp.Variable(shape=preference.shape, boolean=True)
        bind_2 = cp.Variable(shape=preference.shape[1], boolean=True)
        bind_3 = cp.Variable(shape=preference.shape[1], boolean=True)
        bind_constraint = bind_2 + bind_3 == 1

        for vgs in valid_group_sizes:
            min_group_size = min(vgs)
            max_group_size = max(vgs)

            groupmax = np.array([max_group_size] * preference.shape[1])

            group_constraint_1 = cp.sum(selection, axis=0) <= groupmax
            group_constraint_2 = (1 - bind_2) * min_group_size >= min_group_size - cp.sum(selection, axis=0)
            group_constraint_3 = (1 - bind_3) * max_group_size >= cp.sum(selection, axis=0)
            assignment_constraint = cp.sum(selection, axis=1) == 1

            cost = cp.sum(cp.multiply(preference, selection))

            constraints = [group_constraint_1, group_constraint_2,
                           group_constraint_3, bind_constraint,
                           assignment_constraint]

            assign_prob = cp.Problem(cp.Minimize(cost), constraints)

            assign_prob.solve(solver=cp.GLPK_MI)

            solution = pd.DataFrame(selection.value,
                                    index=preference.index,
                                    columns=preference.columns)

            groups = {}
            for index, row in solution.T.iterrows():
                students = list(row[row == 1].index)
                if len(students) > 0:
                    groups[index] = students

            solutions.append(groups)
        return solutions

    @staticmethod
    def sorted_k_partitions(seq, k):
        n = len(seq)
        groups = []  # a list of lists, currently empty

        def generate_partitions(i):
            if i >= n:
                yield list(map(tuple, groups))
            else:
                if n - i > k - len(groups):
                    for group in groups:
                        group.append(seq[i])
                        yield from generate_partitions(i + 1)
                        group.pop()

                if len(groups) < k:
                    groups.append([seq[i]])
                    yield from generate_partitions(i + 1)
                    groups.pop()

        result = generate_partitions(0)

        return result

    @staticmethod
    def score_group(groups, data):
        info = data[0]
        preference = data[1]

        scores = []
        diversity_score = 0
        for group in groups:
            diversity_score += len(info.loc[list(group), "Bachelor"].unique())

        # try each project for each group.
        group_project_combi = [list(zip(x, groups)) for x in itertools.permutations(preference.columns, len(groups))]

        for option in group_project_combi:
            preference_score = 0
            group_assignment = []
            for project, group in option:
                preference_score += preference.loc[list(group), project].sum()
                group_assignment.append((project, group))

            # Save.
            scores.append([group_assignment, preference_score, diversity_score])

        return scores

    def print_arguments(self):
        print("Arguments:")
        print("  > Data file: {}".format(self.data_path))
        print("  > Sheet name: {}".format(self.sheet_name))
        print("  > Min group size: {}".format(self.min_group_size))
        print("  > Max group size: {}".format(self.max_group_size))
        print("  > Cores: {}".format(self.cores))
        print("  > Max. end time: {}".format(datetime.fromtimestamp(self.max_end_time).strftime("%Y-%m-%d %H:%M:%S")))
        print("")


class MultiProcessor:
    def __init__(self, cores, max_end_time, data):
        self.cores = cores
        self.max_end_time = max_end_time
        self.data = data

    def process(self, work_func, iterable):
        cores = self.cores
        if cores > len(iterable):
            cores = len(iterable)

        if cores == 1:
            result = self.singe_process(work_func=work_func,
                                        iterable=iterable)
        else:
            result = self.multiprocessing(work_func=work_func,
                                          iterable=iterable,
                                          cores=cores)

        return result

    def singe_process(self, work_func, iterable):
        print("\t\tProcessing {} instances with 1 core.".format(len(iterable)))

        result = []
        best_preference_score = np.inf
        start_time = int(time.time())
        last_print_time = int(time.time())
        for i, instance in enumerate(iterable):
            if self.max_end_time <= time.time():
                break

            output = work_func(instance, self.data.copy())
            if isinstance(output, list):
                for out in output:
                    if out[1] <= (best_preference_score + 1):
                        best_preference_score = out[1]
                        result.append(out)
            else:
                if output[1] <= (best_preference_score + 1):
                    best_preference_score = output[1]
                    result.append(output)

            now_time = int(time.time())
            if now_time - last_print_time >= 60:
                last_print_time = now_time
                rt_min, rt_sec = divmod(time.time() - start_time, 60)
                rt_hour, rt_min = divmod(rt_min, 60)
                print("\t\t\t[{:02d}:{:02d}:{:02d}] {}/{} instances completed [{:.2f}%]".format(int(rt_hour), int(rt_min), int(rt_sec), i+1, len(iterable), (100 / len(iterable)) * (i+1)))

        return result

    def multiprocessing(self, work_func, iterable, cores):
        # Create queues.
        in_manager = mp.Manager()
        input_q = in_manager.Queue(maxsize=len(iterable))
        out_manager = mp.Manager()
        output_q = out_manager.Queue(maxsize=len(iterable))

        # Populate queues.
        for instance in iterable:
            input_q.put(instance)

        # Create workers.
        processes = []
        for proc_id in range(cores):
            processes.append(mp.Process(target=self.worker,
                                        args=(input_q, output_q, work_func,
                                              self.data.copy(),
                                              self.max_end_time)
                                        ))

        # Start workers.
        started_procs = []
        for proc in processes:
            try:
                proc.start()
                started_procs.append(proc)
            except RuntimeError:
                pass
        del processes

        # Log the progress
        print("\t\tProcessing {} instances with {} cores.".format(len(iterable), cores))
        n = 0
        start_time = int(time.time())
        last_print_time = int(time.time())
        while n < len(iterable):
            time.sleep(2)
            n = output_q.qsize()
            if self.max_end_time <= time.time():
                print("\t\t\tmax progress time reached")
                break

            now_time = int(time.time())
            if now_time - last_print_time >= 60:
                last_print_time = now_time
                rt_min, rt_sec = divmod(time.time() - start_time, 60)
                rt_hour, rt_min = divmod(rt_min, 60)
                print("\t\t\t[{:02d}:{:02d}:{:02d}] {}/{} instances completed [{:.2f}%]".format(int(rt_hour), int(rt_min), int(rt_sec), n, len(iterable), (100 / len(iterable)) * n))

        # Stop the cores.
        for _ in range(cores):
            input_q.put("DONE")

        # Join the workers.
        for proc in started_procs:
            proc.join()

        # Join the results.
        result = []
        best_preference_score = np.inf
        while output_q.qsize() > 0:
            (instance, output) = output_q.get(True, timeout=1)
            if isinstance(output, list):
                for out in output:
                    if out[1] <= (best_preference_score + 1):
                        best_preference_score = out[1]
                        result.append(out)
            else:
                if output[1] <= (best_preference_score + 1):
                    best_preference_score = output[1]
                    result.append(output)

        del in_manager, input_q, out_manager, output_q

        return result

    @staticmethod
    def worker(input_queue, output_queue, work_func, data, max_end_time):
        while True:
            if max_end_time <= time.time():
                break
            try:
                if not input_queue.empty():
                    instance = input_queue.get(True, timeout=1)
                    if instance == "DONE":
                        break

                    output_queue.put((instance, work_func(instance, data)))
            except queue.Empty:
                pass

            time.sleep(random.uniform(0.5, 1))


if __name__ == '__main__':
    m = main()
    m.start()
