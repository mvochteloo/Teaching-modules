import pandas as pd
import ast

groups = pd.read_excel('good_group_distributions.xlsx', index_col=0, header=0)

df = pd.read_excel('project_preferences.xlsx', header=0, sheet_name='Sheet1')
df.set_index("Student", inplace=True)
info = df[["Bachelor"]].copy()
preference = df.iloc[:, 1:].copy()

student_to_study_dict = dict(zip(info.index, info.loc[:, "Bachelor"]))

s = []

for i, (_, (group_assignment, preference_score, diversity_score)) in enumerate(groups.iterrows()):
    match = False

    output = "Solution {}:\n".format(i)
    group_assignment = ast.literal_eval(group_assignment)
    group_assignment.sort(key=lambda x: x[0])

    choice_counts = {}

    for project, members in group_assignment:
        # if match is False and s[0] in members and s[1] in members and s[2] in members:
        #     math = True
        if len(group_assignment) == 4:
            match = True

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
        output += "\t  {}: {} {}\n".format(project, ", ".join(text), studies)

    choice_info = []
    for j in range(1, 100):
        if j in choice_counts:
            choice_info.append("{}x{}".format(j, choice_counts[j]))

    output += "\t  Preference score: {}\t\tDiversity score: {}\t\tN Groups: {}\t\tChoice counts: {}".format(preference_score, diversity_score, len(group_assignment), ", ".join(choice_info))

    if match:
        print(output)

    i += 1