#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
using namespace std;

const uint64_t base = 131;

struct Node {
    int left, right; // 1-indexed positions in seq1 (or seq2 for second move)
    char dir;        // '+' or '-' indicating the direction
    Node() : left(0), right(0), dir(' ') {}
    Node(int l, int r, char d) : left(l), right(r), dir(d) {}
};

struct SubstringEntry {
    uint64_t hash;
    int end; // ending position in seq1 (1-indexed)
};

bool compareSubstringEntry(const SubstringEntry &a, const SubstringEntry &b) {
    if (a.hash != b.hash)
        return a.hash < b.hash;
    return a.end < b.end;
}

vector<uint64_t> computePowers(int n) {
    vector<uint64_t> powers(n + 1, 0);
    powers[0] = 1;
    for (int i = 1; i <= n; i++) {
        powers[i] = powers[i - 1] * base;
    }
    return powers;
}

vector<uint64_t> computePrefixHash(const string &s, const vector<uint64_t> &powers) {
    int n = s.size();
    vector<uint64_t> h(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        h[i] = h[i - 1] * base + static_cast<uint64_t>(s[i - 1]);
    }
    return h;
}

// l, r 1-indexed
uint64_t substringHash(const vector<uint64_t> &prefix, const vector<uint64_t> &powers, int l, int r) {
    return prefix[r] - prefix[l - 1] * powers[r - l + 1];
}

string reverseComplement(const string &s) {
    int n = s.size();
    string rc(n, ' ');
    for (int i = 0; i < n; i++) {
        char comp;
        switch (s[i]) {
            case 'A': comp = 'T'; break;
            case 'T': comp = 'A'; break;
            case 'C': comp = 'G'; break;
            case 'G': comp = 'C'; break;
            default: comp = s[i];
        }
        rc[n - 1 - i] = comp;
    }
    return rc;
}
struct printans{
    int rpos, len, count;
    char dir;
};
printans make_print(int a, int b, int c, char ch)
{
    printans ans;
    ans.rpos = a, ans.len = b, ans.count = c;
    ans.dir = ch;
    return ans;
}
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // input
    string seq1, seq2;
    cin >> seq1 >> seq2;
    int len1 = seq1.size(), len2 = seq2.size();
    // cout << len1 << " " << len2 << endl;
    // the revseq of seq2
    string revComp = reverseComplement(seq2);

    // calc the maximum length
    int maxLen = max(len1, len2);
    vector<uint64_t> powers = computePowers(maxLen);

    // calc hash
    vector<uint64_t> hashSeq1 = computePrefixHash(seq1, powers);
    vector<uint64_t> hashSeq2 = computePrefixHash(seq2, powers);
    vector<uint64_t> hashRevComp = computePrefixHash(revComp, powers);

    vector<SubstringEntry> substrings;
    for (int i = 1; i <= len1; i++) {
        for (int j = i; j <= len1; j++) {
            uint64_t h = substringHash(hashSeq1, powers, i, j);
            substrings.push_back({h, j});
        }
    }
    sort(substrings.begin(), substrings.end(), compareSubstringEntry);
    const int delta = 35;
    // construct rankmap
    unordered_map<uint64_t, int> rankMap;
    if (!substrings.empty()) {
        rankMap[substrings[0].hash] = 1;
    }
    for (size_t i = 1; i < substrings.size(); i++) {
        if (substrings[i].hash != substrings[i - 1].hash) {
            rankMap[substrings[i].hash] = rankMap[substrings[i - 1].hash] + 1;
        }
    }

    // Group the ending positions by rank
    unordered_map<int, vector<int>> entriesByRank;
    for (const auto &entry : substrings) {
        int rank = rankMap[entry.hash];
        entriesByRank[rank].push_back(entry.end);
    }

    // dp[i][j] represents the matchability of seq2 first i chars and seq1 first j chars 
    vector<vector<bool>> dp(len2 + 1, vector<bool>(len1 + 1, false));
    vector<vector<Node>> renew(len2 + 1, vector<Node>(len1 + 1));
    dp[0][0] = true;

    for (int i = 1; i <= len2; i++) {
        // Direct method
        for (int j = 1; j <= len1; j++) {
            if (seq2[i - 1] == seq1[j - 1] && dp[i - 1][j - 1]) {
                dp[i][j] = true;
            }
        }
        // Hash method
        for (int j_idx = i; j_idx >= 1; j_idx--) {
            int segLen = i - j_idx + 1;

            // Case 1: the postive case
            uint64_t currHash = substringHash(hashSeq2, powers, j_idx, i);
            if (rankMap.find(currHash) != rankMap.end()) {
                int rank = rankMap[currHash];
                if (entriesByRank.find(rank) != entriesByRank.end()) {
                    for (int pos : entriesByRank[rank]) {
                        if (pos <= len1 && dp[i - segLen][pos]) {
                            dp[i][pos] = true;
                            renew[i][pos] = Node(pos - segLen + 1, pos, '+');
                        }
                    }
                }
            }

            // Case 2: the negative case
            // Correspondence: The substring [j_idx, i] of seq2 corresponds to the substring [len2 - i + 1, len2 - j_idx + 1] of revComp
            int l_rev = len2 - i + 1, r_rev = len2 - j_idx + 1;
            currHash = substringHash(hashRevComp, powers, l_rev, r_rev);
            if (rankMap.find(currHash) != rankMap.end()) {
                int rank = rankMap[currHash];
                if (entriesByRank.find(rank) != entriesByRank.end()) {
                    for (int pos : entriesByRank[rank]) {
                        if (pos <= len1 && dp[i - segLen][pos]) {
                            dp[i][pos] = true;
                            renew[i][pos] = Node(pos - segLen + 1, pos, '-');
                        }
                    }
                }
            }
        }
    }
    int row = len2, col = len1 - delta;
    vector<Node> moves;
    while (row != 0) {
        if (renew[row][col].left == 0) {
            row--;
            col--;
            continue;
        }
        moves.push_back(renew[row][col]);
        int segLen = renew[row][col].right - renew[row][col].left + 1;
        moves.push_back(Node(row - segLen + 1, row, renew[row][col].dir));
        row -= segLen;
    }

    // print the answer
    vector<printans> Matchseq;
    int refpos = -1, nowlen = 0;
    Matchseq.push_back(make_print(-1, 0, 0, '+'));
    for (int i = moves.size() - 1; i >= 0; i -= 2) {
        int nextlen = moves[i - 1].right - moves[i - 1].left + 1;
        int nextpos = moves[i - 1].right;
        if(refpos != nextpos || nextlen != nowlen){
            Matchseq.push_back(make_print(nextpos, nextlen, 1, moves[i - 1].dir));
            refpos = nextpos, nowlen = nextlen;
        }else{
            Matchseq[Matchseq.size() - 1].count++;
        }
    }
    for(int i = 1; i <= Matchseq.size() - 1; ++ i)
    {
        cout << i << " " << Matchseq[i].rpos << " " << Matchseq[i].len << " " << Matchseq[i].count << " ";
        cout << (Matchseq[i].dir == '+' ? "postive" : "negative") << endl;
    }
    return 0;
}
