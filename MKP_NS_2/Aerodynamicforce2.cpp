#include "Aerodynamicforce2.h"

#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <utility>

// ------------------------------------------------------------
// P2 tvarové funkce na referenčním trojúhelníku
// Vrcholy: A(0,0), B(1,0), C(0,1)
// Barycentrické: L1=1-ξ-η, L2=ξ, L3=η
// Pořadí uzlů (Gmsh 2nd order triangle):
// 0:A, 1:B, 2:C, 3:AB, 4:BC, 5:CA
// ------------------------------------------------------------
static inline void p2_shape_derivs_ref(double xi, double eta, double dN_dxi[6], double dN_deta[6]) {
    const double L1 = 1.0 - xi - eta;
    const double L2 = xi;
    const double L3 = eta;

    // dL/dxi, dL/deta
    const double dL1_dxi = -1.0, dL1_deta = -1.0;
    const double dL2_dxi = 1.0, dL2_deta = 0.0;
    const double dL3_dxi = 0.0, dL3_deta = 1.0;

    // N1 = L1(2L1-1) = 2L1^2 - L1
    // dN1 = (4L1-1) dL1
    const double c1 = (4.0 * L1 - 1.0);
    dN_dxi[0] = c1 * dL1_dxi;
    dN_deta[0] = c1 * dL1_deta;

    // N2 = L2(2L2-1)
    const double c2 = (4.0 * L2 - 1.0);
    dN_dxi[1] = c2 * dL2_dxi;
    dN_deta[1] = c2 * dL2_deta;

    // N3 = L3(2L3-1)
    const double c3 = (4.0 * L3 - 1.0);
    dN_dxi[2] = c3 * dL3_dxi;
    dN_deta[2] = c3 * dL3_deta;

    // N4 = 4 L1 L2
    dN_dxi[3] = 4.0 * (dL1_dxi * L2 + L1 * dL2_dxi);
    dN_deta[3] = 4.0 * (dL1_deta * L2 + L1 * dL2_deta);

    // N5 = 4 L2 L3
    dN_dxi[4] = 4.0 * (dL2_dxi * L3 + L2 * dL3_dxi);
    dN_deta[4] = 4.0 * (dL2_deta * L3 + L2 * dL3_deta);

    // N6 = 4 L3 L1
    dN_dxi[5] = 4.0 * (dL3_dxi * L1 + L3 * dL1_dxi);
    dN_deta[5] = 4.0 * (dL3_deta * L1 + L3 * dL1_deta);
}

// ------------------------------------------------------------
// Klíč pro unordered_map z neorientované hrany (u,v) -> info
// ------------------------------------------------------------
static inline std::uint64_t edge_key(int a, int b) {
    const std::uint32_t lo = (std::uint32_t)std::min(a, b);
    const std::uint32_t hi = (std::uint32_t)std::max(a, b);
    return (std::uint64_t(hi) << 32) | std::uint64_t(lo);
}

struct AdjInfo {
    int tri;     // index trojúhelníku v meshResult.matrix1
    int edgeId;  // 0=AB, 1=BC, 2=CA (lokálně na referenčním trojúhelníku)
    // globální vrcholy trojúhelníku (A,B,C):
    int gA, gB, gC;
    // globální P2 uzly trojúhelníku v gmsh pořadí: A,B,C,AB,BC,CA
    int gP2[6];
};

// ------------------------------------------------------------
// Hlavní funkce: trakce + přesný ∇u z přilehlého trojúhelníku
// ------------------------------------------------------------
std::vector<double> AerodynamicForce2(
    const std::string& /*meshFile*/,
    MatrixVectorResult& meshResult,
    const int& tag,
    const double& /*t*/,
    const std::vector<double>& xvec
) {
    
    const double Re = 50.0;
    const double mu = 2.0 / 3.0 * 0.2 / Re;

    auto& m = meshResult;
    const int nTri = m.integerResult;
    const int N = (int)m.matrix3.size();

    // --- mapování tlaku (shodné s tvým StiffnessMatrix/MassMatrix)
    // tlakové DOF existují jen na "vertex" uzlech (A,B,C z každého tri)
    std::vector<int> isVertex(N, 0);
    for (int k = 0; k < nTri; ++k) {
        int A = m.matrix1[k][0] - 1;
        int B = m.matrix1[k][1] - 1;
        int C = m.matrix1[k][2] - 1;
        isVertex[A] = isVertex[B] = isVertex[C] = 1;
    }
    int trans_i = 0;
    std::vector<int> index_pressure(N, 0);
    for (int i = 0; i < N; ++i) {
        if (isVertex[i]) index_pressure[i] = trans_i++;
    }

    auto p_at_vertex = [&](int node)->double {
        if (!isVertex[node]) return 0.0;
        const int ip = index_pressure[node];
        return xvec[2 * N + ip];
        };

    // --- zrychlení: mapa neorientovaných hran (vrchol,vrchol) -> přilehlý trojúhelník + typ hrany
    std::unordered_map<std::uint64_t, AdjInfo> edge2tri;
    edge2tri.reserve((size_t)nTri * 3);

    for (int k = 0; k < nTri; ++k) {

        const int gA = m.matrix1[k][0] - 1;
        const int gB = m.matrix1[k][1] - 1;
        const int gC = m.matrix1[k][2] - 1;
        const int gAB = m.matrix1[k][3] - 1;
        const int gBC = m.matrix1[k][4] - 1;
        const int gCA = m.matrix1[k][5] - 1;

        AdjInfo info;
        info.tri = k;
        info.gA = gA; info.gB = gB; info.gC = gC;
        info.gP2[0] = gA; info.gP2[1] = gB; info.gP2[2] = gC;
        info.gP2[3] = gAB; info.gP2[4] = gBC; info.gP2[5] = gCA;

        // AB
        info.edgeId = 0;
        edge2tri.emplace(edge_key(gA, gB), info);

        // BC
        info.edgeId = 1;
        edge2tri.emplace(edge_key(gB, gC), info);

        // CA
        info.edgeId = 2;
        edge2tri.emplace(edge_key(gC, gA), info);
    }

    // --- Gauss na [0,1]
    const double Q[3] = {
        0.5 - std::sqrt(3.0 / 5.0) / 2.0,
        0.5,
        0.5 + std::sqrt(3.0 / 5.0) / 2.0
    };
    const double W[3] = { 5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0 };

    auto ux = [&](int node) { return xvec[node]; };
    auto uy = [&](int node) { return xvec[node + N]; };

    double Fx = 0.0, Fy = 0.0;

    // projdi všechny hrany, vezmi jen ty s daným tagem (válec)
    for (int e = 0; e < (int)m.vector1.size(); ++e) {
        if (m.vector1[e] != tag) continue;

        // gmsh 2nd order line: [D,E,F] = [end1,end2,mid]
        const int gD = m.matrix2[e][0] - 1;
        const int gE = m.matrix2[e][1] - 1;
        const int gF = m.matrix2[e][2] - 1;

        // najdi přilehlý trojúhelník podle koncových vrcholů (gD,gE)
        auto it = edge2tri.find(edge_key(gD, gE));
        if (it == edge2tri.end()) {
            // hrana bez přilehlého trojúhelníku v mapě -> nekonzistentní síť / tag
            continue;
        }
        const AdjInfo adj = it->second;

        // trojúhelník (globální uzly)
        const int gA = adj.gP2[0], gB = adj.gP2[1], gC = adj.gP2[2];
        const int gAB = adj.gP2[3], gBC = adj.gP2[4], gCA = adj.gP2[5];

        // geometrie trojúhelníku pro Jacobian J = [B-A, C-A]
        const auto& XA = m.matrix3[gA];
        const auto& XB = m.matrix3[gB];
        const auto& XC = m.matrix3[gC];

        const double J00 = XB[0] - XA[0];
        const double J01 = XC[0] - XA[0];
        const double J10 = XB[1] - XA[1];
        const double J11 = XC[1] - XA[1];

        const double detJ = J00 * J11 - J01 * J10;
        if (std::abs(detJ) < 1e-14) continue;

        // invJ = J^{-1}
        const double invJ00 = J11 / detJ;
        const double invJ01 = -J01 / detJ;
        const double invJ10 = -J10 / detJ;
        const double invJ11 = J00 / detJ;

        // jednotková normála hrany z (gD -> gE) (orientaci můžeš případně flipnout)
        const auto& XD = m.matrix3[gD];
        const auto& XE = m.matrix3[gE];
        const double ex = XE[0] - XD[0];
        const double ey = XE[1] - XD[1];
        const double L = std::sqrt(ex * ex + ey * ey);
        if (L <= 0) continue;
        const double nx = ey / L;
        const double ny = -ex / L;

        // Zjisti, která lokální hrana to je a jak parametrizovat (s=0 na gD, s=1 na gE)
        // Hrany referenčního trojúhelníku:
        // AB: (xi,eta)=(s,0)     od A do B
        // BC: (xi,eta)=(1-s,s)   od B do C
        // CA: (xi,eta)=(0,1-s)   od C do A  (tj. s=0 na C, s=1 na A)
        //
        // My chceme orientaci podle (gD -> gE), takže případně s nahradíme (1-s).

        auto map_edge_to_ref = [&](double s, double& xi, double& eta) {
            if (adj.edgeId == 0) {
                // AB
                const bool D_is_A = (gD == gA);
                const bool D_is_B = (gD == gB);
                // očekáváme {gD,gE} = {gA,gB} v nějakém pořadí
                const double ss = (D_is_A ? s : (1.0 - s));
                xi = ss; eta = 0.0;
            }
            else if (adj.edgeId == 1) {
                // BC: B->C
                const bool D_is_B = (gD == gB);
                const double ss = (D_is_B ? s : (1.0 - s));
                xi = 1.0 - ss;
                eta = ss;
            }
            else {
                // CA: C->A (na referenci)
                // Pokud gD==gC, pak ss=s. Pokud gD==gA, musíme obrátit (1-s).
                const bool D_is_C = (gD == gC);
                const double ss = (D_is_C ? s : (1.0 - s));
                xi = 0.0;
                eta = 1.0 - ss;
            }
            };

        // tlak na hraně: bereme lineárně mezi koncovými vrcholy (gD,gE) – musí to být vertex uzly
        const double pD = p_at_vertex(gD);
        const double pE = p_at_vertex(gE);

        // rychlosti v P2 uzlech trojúhelníku v pořadí [A,B,C,AB,BC,CA]
        const int gU[6] = { gA, gB, gC, gAB, gBC, gCA };
        double uval[6], vval[6];
        for (int i = 0; i < 6; ++i) {
            uval[i] = ux(gU[i]);
            vval[i] = uy(gU[i]);
        }

        // kvadratura po hraně
        for (int q = 0; q < 3; ++q) {
            const double s = Q[q];
            const double w = W[q];

            // (xi,eta) bod na hraně v referenčním trojúhelníku (orientace dle gD->gE)
            double xi, eta;
            map_edge_to_ref(s, xi, eta);

            // ∂Ni/∂xi, ∂Ni/∂eta
            double dN_dxi[6], dN_deta[6];
            p2_shape_derivs_ref(xi, eta, dN_dxi, dN_deta);

            // ∂Ni/∂x, ∂Ni/∂y : grad = invJ^T * grad_ref
            // dN/dx = invJ00*dN/dxi + invJ10*dN/deta
            // dN/dy = invJ01*dN/dxi + invJ11*dN/deta
            double dux_dx = 0.0, dux_dy = 0.0;
            double duy_dx = 0.0, duy_dy = 0.0;

            for (int i = 0; i < 6; ++i) {
                const double dNdx = invJ00 * dN_dxi[i] + invJ10 * dN_deta[i];
                const double dNdy = invJ01 * dN_dxi[i] + invJ11 * dN_deta[i];

                dux_dx += uval[i] * dNdx;
                dux_dy += uval[i] * dNdy;

                duy_dx += vval[i] * dNdx;
                duy_dy += vval[i] * dNdy;
            }

            // tlak v bodě (lineární mezi gD a gE podle s v orientaci gD->gE)
            const double p = pD * (1.0 - s) + pE * s;

            // 2*mu*D(u)
            const double sxx = 2.0 * mu * dux_dx;
            const double syy = 2.0 * mu * duy_dy;
            const double sxy = mu * (dux_dy + duy_dx);

            // trakce t = -p n + (2mu D) n
            const double tx = -p * nx + sxx * nx + sxy * ny;
            const double ty = -p * ny + sxy * nx + syy * ny;

            Fx += tx * (L * w);
            Fy += ty * (L * w);
        }
    }

    return { Fx, Fy };
}
