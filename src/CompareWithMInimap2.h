#ifndef COMPAREWITHMINIMAP2_H
#define COMPAREWITHMINIMAP2_H

#include <unordered_map>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include "Aligner.h"
#include "Alignment.h"
#include "Utils.h"

#define CHECK_EXACT_MATCH_CRITERION__ 0

using namespace std;

class ComparatorWithMinimap2 {
private:
    double percentageIdentity;

public: 
    ComparatorWithMinimap2(double pi) : percentageIdentity(pi) {}

    string convertCigarToStr(const string& cigar, bool skipClips = false) {
        stringstream simplifiedCigar;
        int len = cigar.length();

        for (int i = 0; i < len; ++i) {
            // Parse the numeric part of CIGAR operation
            int num = 0;
            while (i < len && isdigit(cigar[i])) {
                num = num * 10 + (cigar[i] - '0');
                ++i;
            }

            // Extract the CIGAR operation
            char op = cigar[i];

            // Perform actions based on CIGAR operation
            switch (op) {
                case 'M':
                    // Match or mismatch (substitution)
                    simplifiedCigar << string(num, '=');
                    break;
                case 'X':
                    // Match or mismatch (substitution)
                    simplifiedCigar << string(num, 'X');
                    break;
                case 'I':
                    // Insertion
                    simplifiedCigar << string(num, 'I');
                    break;
                case 'D':
                    // Deletion
                    simplifiedCigar << string(num, 'D');
                    break;
                case 'S':
                case 'H':
                    // Soft/Hard clipping
                    if (!skipClips) simplifiedCigar << string(num, 'S');
                    break;
                default:
                    // Handle other CIGAR operations if needed
                    break;
            }
        }

        return simplifiedCigar.str();
    }

    bool hasConsecutiveMatches(const std::string& s, int k) {
        int count = 0;
        for (char c : s) {
            if (c == 'M') {
                count++;
                if (count >= k) {
                    return true;
                }
            } else {
                count = 0; // Reset count if the current character is not 'M'
            }
        }
        return false;
    }

    uint32_t getMatchesCount(string cigarStr) {
        uint32_t cnt = 0;
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'M') {
                cnt++;
            }
        }
        return cnt;
    }

    string trimClips(string cigarStr) {
        string trimmedCigarStr = "";
        for (size_t i = 0; i < cigarStr.length(); i++) {
            if (cigarStr[i] == 'S' || cigarStr[i] == 'H') {
                continue;
            }
            trimmedCigarStr += cigarStr[i];
        }
        return trimmedCigarStr;
    }

    bool checkIdentityPercentange(const Alignment &aln) {
        double identity = (double) aln.matches / (double) aln.partialMatchSize;
        if (DEBUG_MODE) cout << "matches: " << aln.matches << ", cigarStr.size : " << aln.partialMatchSize << ", identity: " << identity << endl;
        if (identity < percentageIdentity)
            return false;
        return true;
    }

    unordered_map<uint16_t, uint32_t> countOccurrences(const multiset<uint16_t>& ms) {
        unordered_map<uint16_t, uint32_t> countMap;

        for (const uint16_t& element : ms) {
            countMap[element]++;
        }

        return countMap;
    }

    void comparePmWithMinimap2(const Config& cfg, tsl::robin_map <uint32_t, string>& reads, tsl::robin_map <uint32_t, string>& queries, const uint32_t queryCount, const string& outputAddress)
    {
        Utilities<uint16_t> utils;
        vector<Alignment> minimap2Alignments, parmikAlignments;

        // Read the Minimap2 SAM file
        SamReader minimap2Sam(cfg.otherToolOutputFileAddress);
        minimap2Sam.parseFile(queryCount, minimap2Alignments, false, true);

        // Read the Parmik SAM file
        SamReader parmikSam(cfg.baselineBaseAddress);
        parmikSam.parseFile(queryCount, parmikAlignments, false, false);

        ofstream outputFile(outputAddress);
        if (!outputFile.is_open()) {
            cout << "error opening file: " << outputAddress << endl;
            return;
        }

        uint32_t minimap2FN = 0, parmikFN = 0, totalTN = 0, numberOfQueryContainN = 0;
        uint64_t minimap2Outperform = 0, parmikOutperform = 0, equalPerformance = 0;
        uint64_t sameReadMinimap2Outperform = 0, differentReadMinimap2Outperform = 0, sameReadPmOutperform = 0, differentReadPmOutperform = 0, sameReadEqual = 0, differentReadEqual = 0;
        uint64_t minimap2LowPI_FN = 0, minimap2LowPIOutperformBestPm = 0;
        
        multiset<uint16_t> sameReadMinimap2OutperformBps, differentReadMinimap2OutperformBps, sameReadPmOutperformBps, differentReadPmOutperformBps, sameReadminimap2LowPIOutperformBps, differentReadminimap2LowPIOutperformBps;
        unordered_map<uint16_t, uint32_t> pmBestAlnSizeWhenMinimap2FN;

        for (uint32_t queryInd = 0; queryInd < queryCount; queryInd++)
        {
            outputFile << "------------------------------------------------------" << endl;
            vector<Alignment> minimap2LowPI_Alignments;

            if (queryInd % 1000 == 0) {
                cout << "queries processed: " << queryInd << " / " << queryCount << endl;
            }

            string query = queries[queryInd];
            if (query.find('N') != string::npos || query.find('n') != string::npos)
            {
                cout << "query contains N!" << endl;
                numberOfQueryContainN++; 
                continue;
            }

            vector<Alignment> query_minimap2Alignments;
            for (const Alignment& aln : minimap2Alignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    if (checkIdentityPercentange(aln)) {
                        if (CHECK_EXACT_MATCH_CRITERION__) {
                            string cigarStr = convertCigarToStr(aln.cigar, true);
                            if (hasConsecutiveMatches(cigarStr, cfg.kmerLength))
                                query_minimap2Alignments.push_back(aln);
                        } else {
                            query_minimap2Alignments.push_back(aln);
                        }
                    } else {
                        minimap2LowPI_Alignments.push_back(aln);
                        minimap2LowPI_FN++;
                    }
                }
            }
            size_t minimap2ReadPerQuery = query_minimap2Alignments.size();

            vector<Alignment> query_parmikAlignments;
            for (const Alignment& aln : parmikAlignments) 
            {
                if ((uint32_t)aln.queryID == queryInd) {
                    if (CHECK_EXACT_MATCH_CRITERION__) {
                        string cigarStr = convertCigarToStr(aln.cigar, true);
                        if (hasConsecutiveMatches(cigarStr, cfg.kmerLength))
                            query_parmikAlignments.push_back(aln);
                    } else {
                        query_parmikAlignments.push_back(aln);
                    }
                }
            }
            size_t parmikReadPerQuery = query_parmikAlignments.size();

            bool bothHasTP = false, minimap2HasLowPIFN = false;
            vector<uint32_t> minimap2TPReadIds;

            if (minimap2ReadPerQuery == 0 && parmikReadPerQuery == 0) {
                // TN
                totalTN++;
            } else if (minimap2ReadPerQuery == 0 && parmikReadPerQuery > 0) {
                // FN
                minimap2FN++;
                if (minimap2LowPI_Alignments.size() > 0) {
                    minimap2HasLowPIFN = true;
                }

                Alignment bestAlnPm;
                for (auto it = query_parmikAlignments.begin(); it != query_parmikAlignments.end(); it++) {
                    Alignment pmAlnn = (*it);
                    if (pmAlnn.matches + pmAlnn.inDels + pmAlnn.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        bestAlnPm = pmAlnn;
                    } else if (pmAlnn.matches + pmAlnn.inDels + pmAlnn.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        if (pmAlnn.inDels + pmAlnn.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) {
                            bestAlnPm = pmAlnn;
                        }
                    }
                }
                uint32_t pmBestAlnSize = bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions;
                pmBestAlnSizeWhenMinimap2FN[pmBestAlnSize]++;
            } else if (minimap2ReadPerQuery > 0 && parmikReadPerQuery == 0) {
                // FN
                parmikFN++;
            } else { // both > 0
                bothHasTP = true;
            }

            if (bothHasTP || minimap2HasLowPIFN) {
                Alignment bestAlnMinimap2;
                Alignment bestAlnPm;

                if (!minimap2HasLowPIFN) {
                    for (auto it = query_minimap2Alignments.begin(); it != query_minimap2Alignments.end(); it++) {
                        Alignment minimap2Alnn = (*it);
                        if (minimap2Alnn.matches + minimap2Alnn.inDels + minimap2Alnn.substitutions > bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions) {
                            bestAlnMinimap2 = minimap2Alnn;
                        } else if (minimap2Alnn.matches + minimap2Alnn.inDels + minimap2Alnn.substitutions == bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions) {
                            if (minimap2Alnn.inDels + minimap2Alnn.substitutions < bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions) {
                                bestAlnMinimap2 = minimap2Alnn;
                            }
                        }
                    }

                    Alignment parmikAlnForMinimap2SameReadID;
                    for (auto it = query_parmikAlignments.begin(); it != query_parmikAlignments.end(); it++) {
                        Alignment pmAlnn = (*it);
                        if (bestAlnMinimap2.matches > 0 && bestAlnMinimap2.readID == pmAlnn.readID) {
                            parmikAlnForMinimap2SameReadID = pmAlnn;
                        }
                        if (pmAlnn.matches + pmAlnn.inDels + pmAlnn.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                            bestAlnPm = pmAlnn;
                        } else if (pmAlnn.matches + pmAlnn.inDels + pmAlnn.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                            if (pmAlnn.inDels + pmAlnn.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) {
                                bestAlnPm = pmAlnn;
                            }
                        }
                    }

                    bool foundSameRead = (bestAlnMinimap2.readID == bestAlnPm.readID);

                    if (bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        minimap2Outperform++;
                        if (foundSameRead) {
                            sameReadMinimap2Outperform++;
                            sameReadMinimap2OutperformBps.insert(bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                            outputFile << "sameReadMinimap2Outperform";
                        } else {
                            differentReadMinimap2Outperform++;
                            differentReadMinimap2OutperformBps.insert(bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                            if (parmikAlnForMinimap2SameReadID.partialMatchSize > 0) outputFile << "parmikAlnForMinimap2SameReadID: " << parmikAlnForMinimap2SameReadID.cigar << ", M: " << parmikAlnForMinimap2SameReadID.matches << ", S: " << parmikAlnForMinimap2SameReadID.substitutions << ", InDels: " << parmikAlnForMinimap2SameReadID.inDels << endl;
                            outputFile << "differentReadMinimap2Outperform";
                        }
                        outputFile << " (larger aln_length) for [" << bestAlnMinimap2.queryID << ", " << bestAlnMinimap2.readID << "]:, minimap2 cigar: " 
                                   << bestAlnMinimap2.cigar << ", MD: " << bestAlnMinimap2.mismatchPositions << ", M: " << bestAlnMinimap2.matches << ", S: " << bestAlnMinimap2.substitutions << ", InDels: " << bestAlnMinimap2.inDels
                                   << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                    } else if (bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions < bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        parmikOutperform++;
                        if (foundSameRead) {
                            outputFile << "sameReadPmOutperform";
                            sameReadPmOutperformBps.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions - (bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions));
                            sameReadPmOutperform++;
                        } else {
                            if (parmikAlnForMinimap2SameReadID.partialMatchSize > 0) outputFile << "parmikAlnForMinimap2SameReadID: " << parmikAlnForMinimap2SameReadID.cigar << ", M: " << parmikAlnForMinimap2SameReadID.matches << ", S: " << parmikAlnForMinimap2SameReadID.substitutions << ", InDels: " << parmikAlnForMinimap2SameReadID.inDels << endl;
                            outputFile << "differentReadPmOutperform";
                            differentReadPmOutperformBps.insert(bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions - (bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions));
                            differentReadPmOutperform++;
                        }
                        outputFile << " (larger aln_length) for [" << bestAlnMinimap2.queryID << ", " << bestAlnMinimap2.readID << "]:, minimap2 cigar: " 
                                   << bestAlnMinimap2.cigar << ", MD: " << bestAlnMinimap2.mismatchPositions << ", M: " << bestAlnMinimap2.matches << ", S: " << bestAlnMinimap2.substitutions << ", InDels: " << bestAlnMinimap2.inDels
                                   << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                    } else if (bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        if (bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) {
                            minimap2Outperform++;
                            if (foundSameRead) {
                                sameReadMinimap2Outperform++;
                                sameReadMinimap2OutperformBps.insert((bestAlnPm.inDels + bestAlnPm.substitutions) - (bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions));
                                outputFile << "sameReadMinimap2Outperform";
                            } else {
                                differentReadMinimap2Outperform++;
                                differentReadMinimap2OutperformBps.insert((bestAlnPm.inDels + bestAlnPm.substitutions) - (bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions));
                                if (parmikAlnForMinimap2SameReadID.partialMatchSize > 0) outputFile << "parmikAlnForMinimap2SameReadID: " << parmikAlnForMinimap2SameReadID.cigar << ", M: " << parmikAlnForMinimap2SameReadID.matches << ", S: " << parmikAlnForMinimap2SameReadID.substitutions << ", InDels: " << parmikAlnForMinimap2SameReadID.inDels << endl;
                                outputFile << "differentReadMinimap2Outperform";
                            }
                            outputFile << " (fewer edits) for [" << bestAlnMinimap2.queryID << ", " << bestAlnMinimap2.readID << "]:, minimap2 cigar: " 
                                       << bestAlnMinimap2.cigar << ", MD: " << bestAlnMinimap2.mismatchPositions << ", M: " << bestAlnMinimap2.matches << ", S: " << bestAlnMinimap2.substitutions << ", InDels: " << bestAlnMinimap2.inDels
                                       << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                        } else if (bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions > bestAlnPm.inDels + bestAlnPm.substitutions) {
                            parmikOutperform++;
                            if (foundSameRead) {
                                outputFile << "sameReadPmOutperform";
                                sameReadPmOutperformBps.insert((bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions) - (bestAlnPm.inDels + bestAlnPm.substitutions));
                                sameReadPmOutperform++;
                            } else {
                                if (parmikAlnForMinimap2SameReadID.partialMatchSize > 0) outputFile << "parmikAlnForMinimap2SameReadID: " << parmikAlnForMinimap2SameReadID.cigar << ", M: " << parmikAlnForMinimap2SameReadID.matches << ", S: " << parmikAlnForMinimap2SameReadID.substitutions << ", InDels: " << parmikAlnForMinimap2SameReadID.inDels << endl;
                                outputFile << "differentReadPmOutperform";
                                differentReadPmOutperformBps.insert((bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions) - (bestAlnPm.inDels + bestAlnPm.substitutions));
                                differentReadPmOutperform++;
                            }
                            outputFile << " (fewer edits) for [" << bestAlnMinimap2.queryID << ", " << bestAlnMinimap2.readID << "]:, minimap2 cigar: " 
                                       << bestAlnMinimap2.cigar << ", MD: " << bestAlnMinimap2.mismatchPositions << ", M: " << bestAlnMinimap2.matches << ", S: " << bestAlnMinimap2.substitutions << ", InDels: " << bestAlnMinimap2.inDels
                                       << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                        } else {
                            equalPerformance++;
                            if (foundSameRead) {
                                sameReadEqual++;
                                outputFile << "sameReadEqual";
                            } else {
                                if ((parmikAlnForMinimap2SameReadID.partialMatchSize > 0) && (bestAlnMinimap2.matches + bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions == parmikAlnForMinimap2SameReadID.matches + parmikAlnForMinimap2SameReadID.inDels + parmikAlnForMinimap2SameReadID.substitutions) && 
                                    (bestAlnMinimap2.inDels + bestAlnMinimap2.substitutions == parmikAlnForMinimap2SameReadID.inDels + parmikAlnForMinimap2SameReadID.substitutions)) {
                                    sameReadEqual++;
                                    outputFile << "sameReadEqual";
                                } else {
                                    if (parmikAlnForMinimap2SameReadID.partialMatchSize > 0) outputFile << "parmikAlnForMinimap2SameReadID: " << parmikAlnForMinimap2SameReadID.cigar << ", M: " << parmikAlnForMinimap2SameReadID.matches << ", S: " << parmikAlnForMinimap2SameReadID.substitutions << ", InDels: " << parmikAlnForMinimap2SameReadID.inDels << endl;
                                    differentReadEqual++;
                                    outputFile << "differentReadEqual";
                                }
                            }
                            outputFile << " for [" << bestAlnMinimap2.queryID << ", " << bestAlnMinimap2.readID << "]:, minimap2 cigar: " 
                                       << bestAlnMinimap2.cigar << ", MD: " << bestAlnMinimap2.mismatchPositions << ", M: " << bestAlnMinimap2.matches << ", S: " << bestAlnMinimap2.substitutions << ", InDels: " << bestAlnMinimap2.inDels
                                       << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                        }
                    }
                }

                // Check whether Minimap2 Low PI FN alignments are larger than its best alignment
                for (auto it = minimap2LowPI_Alignments.begin(); it != minimap2LowPI_Alignments.end(); it++) {
                    Alignment minimap2LowPIAln = (*it);
                    if (minimap2LowPIAln.matches + minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions > bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        minimap2LowPIOutperformBestPm++;
                        if (minimap2LowPIAln.readID == bestAlnPm.readID) {
                            outputFile << "sameReadminimap2LowPIAlnOutperform";
                            sameReadminimap2LowPIOutperformBps.insert(minimap2LowPIAln.matches + minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                        } else {
                            outputFile << "differentReadminimap2LowPIAlnOutperform";
                            differentReadminimap2LowPIOutperformBps.insert(minimap2LowPIAln.matches + minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions - (bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions));
                        }
                        outputFile << " (larger aln_length) for [" << minimap2LowPIAln.queryID << ", " << minimap2LowPIAln.readID << "]:, minimap2 cigar: " 
                                   << minimap2LowPIAln.cigar << ", MD: " << minimap2LowPIAln.mismatchPositions << ", M: " << minimap2LowPIAln.matches << ", S: " << minimap2LowPIAln.substitutions << ", InDels: " << minimap2LowPIAln.inDels
                                   << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                    } else if (minimap2LowPIAln.matches + minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions == bestAlnPm.matches + bestAlnPm.inDels + bestAlnPm.substitutions) {
                        if (minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions < bestAlnPm.inDels + bestAlnPm.substitutions) {
                            minimap2LowPIOutperformBestPm++;
                            if (minimap2LowPIAln.readID == bestAlnPm.readID) {
                                outputFile << "sameReadminimap2LowPIAlnOutperform";
                                sameReadminimap2LowPIOutperformBps.insert(bestAlnPm.inDels + bestAlnPm.substitutions - (minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions));
                            } else {
                                outputFile << "differentReadminimap2LowPIAlnOutperform";
                                differentReadminimap2LowPIOutperformBps.insert(bestAlnPm.inDels + bestAlnPm.substitutions - (minimap2LowPIAln.inDels + minimap2LowPIAln.substitutions));
                            }
                            outputFile << " (fewer Edits) for [" << minimap2LowPIAln.queryID << ", " << minimap2LowPIAln.readID << "]:, minimap2 cigar: " 
                                       << minimap2LowPIAln.cigar << ", MD: " << minimap2LowPIAln.mismatchPositions << ", M: " << minimap2LowPIAln.matches << ", S: " << minimap2LowPIAln.substitutions << ", InDels: " << minimap2LowPIAln.inDels
                                       << " - parmik cigar: " << bestAlnPm.cigar << ", M: " << bestAlnPm.matches << ", S: " << bestAlnPm.substitutions << ", InDels: " << bestAlnPm.inDels << endl;
                        }
                    }
                }
            }
        }

        outputFile << "<<<<<<<<<<<<<<<<<<<<<<<<Final Results>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        outputFile << "Total Minimap2 TN : " << totalTN << endl;
        outputFile << "Total Minimap2 FN : " << minimap2FN << endl;
        outputFile << "Total Parmik FNN : " << parmikFN << endl;
        outputFile << "Total Minimap2 Outperform : " << minimap2Outperform << endl;
        outputFile << "Total Parmik Outperform : " << parmikOutperform << endl;
        outputFile << "Total Minimap2 Equal : " << equalPerformance << endl;
        outputFile << "Same Read Minimap2 Outperform : " << sameReadMinimap2Outperform << endl;
        outputFile << "Different Read Minimap2 Outperform : " << differentReadMinimap2Outperform << endl;
        outputFile << "Same Read Parmik Outperform : " << sameReadPmOutperform << endl;
        outputFile << "Different Read Parmik Outperform : " << differentReadPmOutperform << endl;
        outputFile << "Same Read Equal : " << sameReadEqual << endl;
        outputFile << "Different Read Equal : " << differentReadEqual << endl;
        outputFile << "minimap2LowPI_FN : " << minimap2LowPI_FN << endl;
        outputFile << "minimap2LowPI_OutperformBestPm : " << minimap2LowPIOutperformBestPm << endl;

        pair<uint16_t, uint16_t> sameReadMinimap2OutperformBpsTuple = utils.calculateStatistics2(sameReadMinimap2OutperformBps);
        outputFile << "No. of Bp (same read) Minimap2 outperforms => [average: " << std::get<0>(sameReadMinimap2OutperformBpsTuple) << ", median: " << std::get<1>(sameReadMinimap2OutperformBpsTuple) << "]" << std::endl;
        
        pair<uint16_t, uint16_t> differentReadMinimap2OutperformBpsTuple = utils.calculateStatistics2(differentReadMinimap2OutperformBps);
        outputFile << "No. of Bp (different read) Minimap2 outperforms => [average: " << std::get<0>(differentReadMinimap2OutperformBpsTuple) << ", median: " << std::get<1>(differentReadMinimap2OutperformBpsTuple) << "]" << std::endl;
        
        pair<uint16_t, uint16_t> sameReadPmOutperformBpsTuple = utils.calculateStatistics2(sameReadPmOutperformBps);
        outputFile << "No. of Bp (same read) PARMIK outperforms => [average: " << std::get<0>(sameReadPmOutperformBpsTuple) << ", median: " << std::get<1>(sameReadPmOutperformBpsTuple) << "]" << std::endl;
        
        pair<uint16_t, uint16_t> differentReadPmOutperformBpsTuple = utils.calculateStatistics2(differentReadPmOutperformBps);
        outputFile << "No. of Bp (different read) PARMIK outperforms => [average: " << std::get<0>(differentReadPmOutperformBpsTuple) << ", median: " << std::get<1>(differentReadPmOutperformBpsTuple) << "]" << std::endl;
        
        pair<uint16_t, uint16_t> sameReadminimap2LowPIOutperformBpsTuple = utils.calculateStatistics2(sameReadminimap2LowPIOutperformBps);
        outputFile << "No. of Bp (same read) Minimap2 outperforms (LowPI) => [average: " << std::get<0>(sameReadminimap2LowPIOutperformBpsTuple) << ", median: " << std::get<1>(sameReadminimap2LowPIOutperformBpsTuple) << "]" << std::endl;
        
        pair<uint16_t, uint16_t> differentReadminimap2LowPIOutperformBpsTuple = utils.calculateStatistics2(differentReadminimap2LowPIOutperformBps);
        outputFile << "No. of Bp (different read) Minimap2 outperforms (LowPI) => [average: " << std::get<0>(differentReadminimap2LowPIOutperformBpsTuple) << ", median: " << std::get<1>(differentReadminimap2LowPIOutperformBpsTuple) << "]" << std::endl;

        // Output histogram of differences
        outputFile << "<<<<<<<<<<<<<Histogram of differences>>>>>>>>>>>>>" << endl;
        outputFile << "Minimap2 > PARMIK (same read)" << endl;
        unordered_map<uint16_t, uint32_t> sameReadMinimap2OutperformMap = countOccurrences(sameReadMinimap2OutperformBps);
        for (const auto& pair : sameReadMinimap2OutperformMap) {
            outputFile << pair.first << ": " << pair.second << endl;
        }

        outputFile << "Minimap2 > PARMIK (different read)" << endl;
        unordered_map<uint16_t, uint32_t> differentReadMinimap2OutperformMap = countOccurrences(differentReadMinimap2OutperformBps);
        for (const auto& pair : differentReadMinimap2OutperformMap) {
            outputFile << pair.first << ": " << pair.second << endl;
        }

        outputFile << "PARMIK > Minimap2 (same read)" << endl;
        unordered_map<uint16_t, uint32_t> sameReadPmOutperformMap = countOccurrences(sameReadPmOutperformBps);
        for (const auto& pair : sameReadPmOutperformMap) {
            outputFile << pair.first << ": " << pair.second << endl;
        }

        outputFile << "PARMIK > Minimap2 (different read)" << endl;
        unordered_map<uint16_t, uint32_t> differentReadPmOutperformMap = countOccurrences(differentReadPmOutperformBps);
        for (const auto& pair : differentReadPmOutperformMap) {
            outputFile << pair.first << ": " << pair.second << endl;
        }

        outputFile << "<<<<<<<<<<<<<<<<<<<<<<<<PARMIK best aln size histo when Minimap2 FN>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        for (const auto& pair : pmBestAlnSizeWhenMinimap2FN) {
            outputFile << pair.first << ": " << pair.second << endl;
        }

        outputFile.close();
    }
};

#endif