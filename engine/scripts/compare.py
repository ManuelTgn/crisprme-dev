
# pos: 232683107, cigarx: XX==X=========I===D=X, seq: AAAAAAAAGTTGCTTTGTATAT
mined_positions = set()
with open("mine_chr1_positive_pos.txt") as fa:
    for i, line in enumerate(fa, 1):
        if "pos:" in line:

            position, cigarx, seq = line.split(",")
            position = position.split("pos:")[1].strip()
            cigarx = cigarx.split("cigarx:")[1].strip()
            seq = seq.split("seq:")[1].strip()

            mined_positions.add(position)


total_correct_positions = 0
total_matched = 0

# Check that we mined at least all the correct off-targets
# RNA,DNA	CTAACAGTTGCTT-TTATCACNNN	CcAtCAGTTGCTTCTT-caACAGG	chr1	231483770	231483770	-	4	2	6
with open("offT_hg38.targets.txt") as fb:
    for line in fb:
        parts = line.strip().split('\t')
        if parts[6] == '+':
            total_correct_positions += 1
            position = parts[5]

            if position in mined_positions:
                print(f"✅ pos {position}")
                total_matched += 1
            else:

                found = False
                deltas = [ -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
                for delta in deltas:
                    position_near = str(int(position) + delta)
                    if not found and position_near in mined_positions:
                        print(f"⚠️ pos {position} (mined {position_near})")
                        total_matched += 1
                        found = True

                if not found:
                    print(f"❌ pos {position}")

percentage_found = total_matched / total_correct_positions * 100
print(f"found {total_matched}/{total_correct_positions} correct positions ({percentage_found}%)")


