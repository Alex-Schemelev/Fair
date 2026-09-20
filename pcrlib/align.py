def simple_align(seq1: str, seq2: str):
    # визуализируем комплементарность: seq1 vs reverse-complement(seq2)
    s1 = seq1.upper()
    # reverse-complement с поддержкой IUPAC
    comp = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B",
        "D": "H", "H": "D", "N": "N"
    }
    rc2 = "".join(comp.get(b, "N") for b in seq2.upper()[::-1])
    s2 = rc2
    best_offset = 0
    best_matches = ""
    best_score = -1

    max_shift = len(s1)
    for offset in range(-max_shift, max_shift + 1):
        matches = []
        score = 0
        for i in range(len(s1)):
            j = i - offset
            if 0 <= j < len(s2):
                if s1[i] == s2[j]:
                    matches.append("|")
                    score += 1
                else:
                    matches.append(" ")
            else:
                matches.append(" ")
        if score > best_score:
            best_score = score
            best_offset = max(0, offset)
            best_matches = "".join(matches)
    return best_offset, best_matches
