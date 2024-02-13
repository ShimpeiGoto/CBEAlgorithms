#pragma once
#include <itensor/all_mps.h>
#include "ShrewdSelection.h"
#include "LocalHamiltonian.h"
#include <chrono>
#include <iomanip>

namespace CBE_DMRG {
        inline double CBE_DMRG(itensor::MPS& psi, const itensor::MPO&H, itensor::Args& args) {
                int Ns = itensor::length(psi);
                itensor::LocalMPO PH(H, args);
                double energy = 0.0;
                auto show_summary = args.getBool("ShowSummary", true);
                double cutoff = args.getReal("Cutoff", 1e-10);
                std::chrono::system_clock::time_point start_dmrg, end_dmrg;

                if (show_summary) {
                        start_dmrg = std::chrono::system_clock::now();
                }
                // Left to Right sweep
                for (int left=1; left < Ns; left++) {
                        psi.position(left);
                        int right = left+1;
                        PH.position(left, psi);
                        itensor::Index siteInd = itensor::siteIndex(psi, left);
                        itensor::Index leftLink = itensor::leftLinkIndex(psi, left);
                        long maxdim = std::min(
                                        itensor::dim(siteInd) * itensor::dim(leftLink),
                                        itensor::dim(siteIndex(psi, right)) * itensor::dim(itensor::rightLinkIndex(psi, right))
                                        );
                        auto [A, Lambda, V] = itensor::svd(psi(left), {leftLink, siteInd});
                        itensor::ITensor B = V * psi(right);
                        itensor::ITensor C = A * Lambda;
                        auto LBlink = itensor::commonIndex(B, Lambda);
                        long rightdim = itensor::dim(rightLinkIndex(psi, left));
                        itensor::ITensor AA;
                        if (rightdim < maxdim) {
                                auto Btr = ShrewdSelection::ShrewdSelectionLeftToRight(left, A, Lambda, B, H, PH, args);
                                auto [Bex, Ind] = itensor::directSum(B, Btr,
                                                                     LBlink, itensor::findIndex(Btr, tinyformat::format("l=%d", left)),
                                                                     {"Tags", tinyformat::format("Link,l=%d", left)});
                                itensor::ITensor Cex = B * itensor::dag(Bex);
                                Cex *= C;
                                bool orthogonal_check = itensor::norm(Cex*Bex - C*B) < cutoff;
                                if (orthogonal_check) {
                                        itensor::ITensor HR = (PH.R() ? PH.R() * Bex : Bex);
                                        HR *= H(right);
                                        HR *= itensor::prime(itensor::dag(Bex));
                                        LocalHamiltonian::SingleSiteH singleH(HR, PH.L(), H(left));
                                        energy = davidson(singleH, Cex);
                                        AA = Cex*Bex;
                                } else { // Expansion might fail
                                        itensor::ITensor HR = (PH.R() ? PH.R() * B : B);
                                        HR *= H(right);
                                        HR *= itensor::prime(itensor::dag(B));
                                        LocalHamiltonian::SingleSiteH singleH(HR, PH.L(), H(left));
                                        energy = davidson(singleH, C);
                                        AA = C*B;
                                }
                        } else {
                                itensor::ITensor HR = (PH.R() ? PH.R() * B : B);
                                HR *= H(right);
                                HR *= itensor::prime(itensor::dag(B));
                                LocalHamiltonian::SingleSiteH singleH(HR, PH.L(), H(left));
                                energy = davidson(singleH, C);
                                AA = C*B;
                        }
                        psi.svdBond(left, AA, itensor::Fromleft, args);
                        psi.normalize();
                }

                // Right to Left sweep
                for (int right=Ns-1; right > 1; right--) {
                        int left = right-1;
                        psi.position(right);
                        PH.position(left, psi);
                        itensor::Index siteInd = itensor::siteIndex(psi, left);
                        itensor::Index leftLink = itensor::leftLinkIndex(psi, left);
                        long maxdim = std::min(
                                        itensor::dim(siteInd) * itensor::dim(leftLink),
                                        itensor::dim(siteIndex(psi, right)) * itensor::dim(itensor::rightLinkIndex(psi, right))
                                        );
                        auto [U, Lambda, B] = itensor::svd(psi(right), {itensor::leftLinkIndex(psi, right)});
                        itensor::ITensor A = psi(left) * U;
                        itensor::ITensor C = Lambda * B;
                        auto ALlink = itensor::commonIndex(A, Lambda);
                        long leftdim = itensor::dim(leftLinkIndex(psi, right));
                        itensor::ITensor AA;
                        if (leftdim < maxdim) {
                                auto Atr = ShrewdSelection::ShrewdSelectionRightToLeft(right, A, Lambda, B, H, PH, args);
                                auto [Aex, Ind] = itensor::directSum(A, Atr,
                                                                     ALlink, itensor::findIndex(Atr, tinyformat::format("l=%d", left)),
                                                                     {"Tags", tinyformat::format("Link,l=%d", left)});
                                itensor::ITensor Cex = A * itensor::dag(Aex);
                                Cex *= C;
                                bool orthogonal_check = itensor::norm(Aex*Cex - A*C) < cutoff;
                                if (orthogonal_check) {
                                        itensor::ITensor HL = (PH.L() ? PH.L() * Aex : Aex);
                                        HL *= H(left);
                                        HL *= itensor::prime(itensor::dag(Aex));
                                        LocalHamiltonian::SingleSiteH singleH(HL, PH.R(), H(left+1));
                                        energy = davidson(singleH, Cex);
                                        AA = Aex * Cex;
                                } else { // Expansion might fail
                                        itensor::ITensor HL = (PH.L() ? PH.L() * A : A);
                                        HL *= H(left);
                                        HL *= itensor::prime(itensor::dag(A));
                                        LocalHamiltonian::SingleSiteH singleH(HL, PH.R(), H(left+1));
                                        energy = davidson(singleH, C);
                                        AA = A * C;
                                }
                        } else {
                                itensor::ITensor HL = (PH.L() ? PH.L() * A : A);
                                HL *= H(left);
                                HL *= itensor::prime(itensor::dag(A));
                                LocalHamiltonian::SingleSiteH singleH(HL, PH.R(), H(right));
                                energy = davidson(singleH, C);
                                AA = A * C;
                        }
                        psi.svdBond(left, AA, itensor::Fromright, args);
                        psi.normalize();
                }

                energy = itensor::innerC(psi, H, psi).real();
                if (show_summary) {
                        end_dmrg = std::chrono::system_clock::now();
                        std::streamsize ss = std::cout.precision();
                        std::cout << "Elapsed time:" << std::fixed << std::setprecision(3) << std::chrono::duration_cast<std::chrono::milliseconds>(end_dmrg - start_dmrg).count() / 1000.0 << "s" << std::endl;
                        std::cout << "Obtained energy:" << std::setprecision(16) << energy << std::endl;
                        std::cout << std::setprecision(ss);
                        std::cout << std::defaultfloat;
                }

                return energy;
        }
}
