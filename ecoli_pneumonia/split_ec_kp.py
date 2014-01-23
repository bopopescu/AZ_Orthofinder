from collections import defaultdict

ec_groups = []
kp_groups = []

with open('../orthogroups.tsv') as f:
    ec_group = []
    kp_group = []

    for line in f:
        line = line.strip()
        if not line:
            ec_groups.append(ec_group)
            kp_groups.append(kp_group)
            ec_group = []
            kp_group = []
            continue

        tokens = line.split('\t')
        group_num = int(tokens[0])
        strain = tokens[1]

        if strain.startswith('Escherichia coli'):
            ec_group.append(tokens)

        elif strain.startswith('Klebsiella'):
            kp_group.append(tokens)

        else:
            print 'Warning!', strain


with open('../ecoli.tsv', 'w') as ec:
    print(len(ec_groups))
    i = 1
    for group in ec_groups:
        strains = set(tokens[1] for tokens in group)
        if len(strains) > 1:
            for tokens in group:
                ec.write('\t'.join([str(i)] + tokens[1:]))
                ec.write('\n')
            ec.write('\n')
            i += 1


with open('../kpneumo.tsv', 'w') as kp:
    print(len(kp_groups))
    i = 1
    for group in kp_groups:
        strains = set(tokens[1] for tokens in group)
        if len(strains) > 1:
            for tokens in group:
                kp.write('\t'.join([str(i)] + tokens[1:]))
                kp.write('\n')
            kp.write('\n')
            i += 1



