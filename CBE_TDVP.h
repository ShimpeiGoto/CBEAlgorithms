#pragma once
#include <itensor/all_mps.h>
#include "ShrewdSelection.h"
#include "LocalHamiltonian.h"

namespace CBE_TDVP {
        template<typename scalar>inline double CBE_TDVP_impl(itensor::MPS& psi, const itensor::MPO& H, scalar dt, itensor::Args& args) {
                int Ns = itensor::length(psi);
                itensor::LocalMPO PH(H, args);
                double cutoff = args.getReal("Cutoff", 1e-10);

                // Left to Right sweep
                for (int left=1; left < Ns; left++) {
                        int right = left+1;
                        psi.position(left);
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
                        auto LBlink = itensor::commonIndex(V, Lambda);
                        long rightdim = itensor::dim(rightLinkIndex(psi, left));
                        itensor::ITensor Bex, Cex;
                        bool expand = false;
                        if (rightdim < maxdim) {
                                auto Btr = ShrewdSelection::ShrewdSelectionLeftToRight(left, A, Lambda, B, H, PH, args);
                                auto [Btmp, Ind] = itensor::directSum(B, Btr,
                                                                     LBlink, itensor::findIndex(Btr, tinyformat::format("l=%d", left)),
                                                                     {"Tags", tinyformat::format("Link,l=%d", left)});
                                Cex = B * itensor::dag(Btmp);
                                Cex *= C;
                                expand = itensor::norm(Cex*Btmp - C*B) < cutoff;
                                if (expand) {
                                        Bex = Btmp;
                                }
                        }

                        if (!expand){
                                Cex = C;
                                Bex = B;
                        }

                        itensor::ITensor HR = (PH.R() ? PH.R() * Bex : Bex);
                        HR *= H(right);
                        HR *= itensor::prime(itensor::dag(Bex));
                        LocalHamiltonian::SingleSiteH singleH(HR, PH.L(), H(left));
                        itensor::applyExp(singleH, Cex, 0.5*dt, args);
                        auto lefttags = tinyformat::format("Link,l=%d", left);
                        auto [Anext, Lnext, Bnext] = itensor::svd(Cex, {leftLink, siteInd}, args+itensor::Args({"LeftTags", lefttags}));
                        auto LambdaNext = Lnext * Bnext;
                        psi.ref(left) = Anext;

                        itensor::ITensor HL = (PH.L() ? PH.L() * Anext : Anext);
                        HL *= H(left);
                        HL *= itensor::prime(itensor::dag(Anext));
                        LocalHamiltonian::ZeroSiteH zeroH(HL, HR);
                        itensor::applyExp(zeroH, LambdaNext, -0.5*dt, args);
                        psi.ref(right) = LambdaNext * Bex;
                }

                // Right edge
                {
                        PH.numCenter(1);
                        psi.position(Ns);
                        PH.position(Ns, psi);

                        LocalHamiltonian::SingleSiteH singleH(PH.L(), PH.R(), H(Ns));
                        itensor::ITensor phi = psi(Ns);
                        itensor::applyExp(singleH, phi, dt, args);
                        auto righttags = tinyformat::format("Link,l=%d", Ns-1);
                        auto [A, L, B] = itensor::svd(phi, {itensor::leftLinkIndex(psi, Ns)}, args+itensor::Args({"RightTags", righttags}));
                        psi.ref(Ns) = B;
                        itensor::ITensor LambdaNext = A*L;
                        itensor::ITensor HR = B;
                        HR *= H(Ns);
                        HR *= itensor::prime(itensor::dag(B));
                        LocalHamiltonian::ZeroSiteH zeroH(PH.L(), HR);
                        itensor::applyExp(zeroH, LambdaNext, -0.5*dt, args);
                        psi.ref(Ns-1) *= LambdaNext;
                }

                PH.numCenter(2);

                // Right to Left sweep
                for (int right = Ns-1; right > 1; right--) {
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
                        bool expand = false;
                        itensor::ITensor Aex, Cex;
                        if (leftdim < maxdim) {
                                auto Atr = ShrewdSelection::ShrewdSelectionRightToLeft(right, A, Lambda, B, H, PH, args);
                                auto [Atmp, Ind] = itensor::directSum(A, Atr,
                                                                     ALlink, itensor::findIndex(Atr, tinyformat::format("l=%d", left)),
                                                                     {"Tags", tinyformat::format("Link,l=%d", left)});
                                Cex = A * itensor::dag(Atmp);
                                Cex *= C;
                                expand = itensor::norm(Atmp*Cex - A*C) < cutoff;
                                if (expand) {
                                        Aex = Atmp;
                                }
                        }

                        if (!expand){
                                Aex = A;
                                Cex = C;
                        }


                        itensor::ITensor HL = (PH.L() ? PH.L() * Aex : Aex);
                        HL *= H(left);
                        HL *= itensor::prime(itensor::dag(Aex));
                        LocalHamiltonian::SingleSiteH singleH(HL, PH.R(), H(right));
                        itensor::applyExp(singleH, Cex, 0.5*dt, args);
                        auto righttags = tinyformat::format("Link,l=%d", left);
                        auto [Anext, Lnext, Bnext] = itensor::svd(Cex, {itensor::commonIndex(Aex, Cex)}, args+itensor::Args({"RightTags", righttags}));
                        auto LambdaNext = Anext * Lnext;
                        psi.ref(right) = Bnext;

                        itensor::ITensor HR = (PH.R() ? PH.R() * Bnext : Bnext);
                        HR *= H(right);
                        HR *= itensor::prime(itensor::dag(Bnext));
                        LocalHamiltonian::ZeroSiteH zeroH(HL, HR);
                        itensor::applyExp(zeroH, LambdaNext, -0.5*dt, args);
                        psi.ref(left) = Aex * LambdaNext;
                }

                // Left edge
                {
                        PH.numCenter(1);
                        psi.position(1);
                        PH.position(1, psi);

                        LocalHamiltonian::SingleSiteH singleH(PH.L(), PH.R(), H(1));
                        itensor::applyExp(singleH, psi.ref(1), 0.5*dt, args);
                }

                return psi.normalize();
        }

        inline double CBE_TDVP(itensor::MPS& psi, const itensor::MPO& H, double dt, itensor::Args& args) {
                std::complex<double> tau(0.0, -dt);
                return CBE_TDVP_impl(psi, H, tau, args);
        }
}
