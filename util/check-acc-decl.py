#!/usr/bin/python3
# Compare variable declaration in 'params.f90' with clauses
# '!$ACC declare' in 'params.f90' and '!$ACC update' in 'setflac.f90'
# Usage: Run this script under src/ directory.

decli, decld, accdecl, accupd = [], [], [], []
acc_hdr = '!$ACC '

s = open('params.f90').readlines()
nline = len(s)
n = 0
while n < nline:
    line = s[n].strip()

    # declaration for integer
    head = 'integer :: '
    if line.startswith(head):
        line = line[len(head):].strip()
        while line.endswith('&'):
            decli.append(line)
            n += 1
            line = s[n].strip()
        decli.append(line) # last line of integer declaration
        #print(line, decli)

    # declaration for double
    head = 'real*8 :: '
    if line.startswith(head):
        line = line[len(head):].strip()
        while line.endswith('&'):
            decld.append(line)
            n += 1
            line = s[n].strip()
        decld.append(line) # last line of double declaration
        #print(line, decld)


    # acc declaration for integer and double
    head = acc_hdr +'declare create('
    if line.startswith(head):
        line = line[len(head):].strip()
        while line.endswith('&'): # multi-lines of acc declare
            accdecl.append(line)
            n += 1
            line = s[n][len(acc_hdr):].strip()
        accdecl.append(line.strip(')')) # last line of acc declare
        #print(line, accdecl)

    n += 1


s = open('setflac.f90').readlines()
nline = len(s)
n = 0
while n < nline:
    line = s[n].strip()

    # acc update for integer and double
    head = acc_hdr +'update device('
    if line.startswith(head):
        line = line[len(head):].strip()
        if not line.endswith('&'): # single line of acc, skip it
            n += 1
            line = s[n].strip()
            continue
        while line.endswith('&'): # multi-lines of acc update
            accupd.append(line)
            n += 1
            line = s[n][len(acc_hdr):].strip()
        accupd.append(line[:-len(') async(1)')])  # last line of acc update
    n += 1
#print(accupd)


# diff'ing
decl = decli + decld
print('Comparing the declaration of F90, ACC-declare, and ACC-update:')
print('# of lines:', len(decl), len(accdecl), len(accupd))

ndiff = 0
for n, line in enumerate(decl):
    if line != accdecl[n]:
        ndiff += 1
        print('Line', n, 'mismatch in acc create:', line, accdecl[n])

    if line != accupd[n]:
        ndiff += 1
        print('Line', n, 'mismatch in acc update:', line, accupd[n])

print('Comparison finished,', ndiff, 'mismatch(es) found.')