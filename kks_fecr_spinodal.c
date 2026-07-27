/* =====================================================================
 * kks_fecr_spinodal.c
 *
 * CALPHAD-informed Kim-Kim-Suzuki (KKS) phase-field model of spinodal
 * decomposition in a binary alloy with a miscibility gap, using
 * thermodynamic data parsed directly from an openly available CALPHAD
 * .tdb file, and a semi-implicit Fourier-spectral (FFTW) time integrator.
 *
 * SYSTEM: Fe-Cr, BCC_A2 (alpha) phase, which has a low-temperature
 *         alpha/alpha' miscibility gap (the classic "885 F embrittlement"
 *         decomposition in ferritic stainless/reactor steels).
 *
 * THERMODYNAMICS (CALPHAD / TDB):
 *   The BCC_A2 molar Gibbs energy is built from parameters parsed at
 *   runtime out of fecr_bcc.tdb:
 *     G = xFe*GHSERFE(T) + xCr*GHSERCR(T)                 (references)
 *       + RT[xCr ln xCr + xFe ln xFe]                     (ideal mixing)
 *       + xCr*xFe*L0(T)                                   (Redlich-Kister
 *                                                            excess)
 *       + RT ln(beta+1) f(T/Tc)                           (Hillert-Jarl /
 *                                                            Inden magnetic
 *                                                            ordering)
 *   fecr_bcc.tdb is a trimmed excerpt of the openly available
 *   multi-component database "mpea-02b.tdb" (B. Hallstedt, 2016/17,
 *   https://github.com/pycalphad/binder), which for the Fe-Cr binary
 *   reproduces the classic assessment of
 *     J.-O. Andersson, B. Sundman, CALPHAD 11 (1987) 83-92.
 *   A tiny hand-written TDB parser (see tdb_parse.h below, inlined here)
 *   reads the FUNCTION/PAR statements and evaluates them at any (x,T) --
 *   nothing about GHSERFE, GHSERCR, or the interaction parameters is
 *   hard-coded as physics; only the *TDB syntax* and the standard CALPHAD
 *   model equations (ideal solution + Redlich-Kister + Hillert-Jarl) are.
 *
 * PHASE-FIELD MODEL -- two selectable modes:
 *
 *  MODE "kks" (default): Kim-Kim-Suzuki, two fields.
 *   Fields: c(r,t)  -- conserved overall Cr mole fraction
 *           eta(r,t) -- non-conserved phase indicator (0 = "alpha" copy,
 *                       1 = "alpha'" copy of the SAME BCC_A2 CALPHAD phase)
 *   KKS constraint at every point:
 *       c = h(eta) cB + (1-h(eta)) cA
 *       dG/dx(cA) = dG/dx(cB)      (equal diffusion potential)
 *   solved by 2x2 Newton iteration using the CALPHAD G(x,T) above.
 *
 *   Free energy functional:
 *     F = INT [ h(eta) G(cB,T) + (1-h(eta)) G(cA,T) + w g(eta)
 *               + (kappa_c/2)|grad c|^2 + (kappa_eta/2)|grad eta|^2 ] dV
 *   h(eta) = eta^3 (6eta^2 -15eta +10),  g(eta) = eta^2 (1-eta)^2
 *
 *   Evolution:
 *     dc/dt   =  div( M grad(mu) ),      mu = dG/dx(cA) [= dG/dx(cB)]
 *     deta/dt = -L [ h'(eta) ( G(cB)-G(cA) - (cB-cA) mu ) + w g'(eta)
 *                    - kappa_eta lap(eta) ]
 *
 *   eta is here just a bookkeeping field (both "phases" share one CALPHAD
 *   description), so it should track c essentially instantaneously. To
 *   guarantee that separation of timescales we set L = R_LM * M with
 *   R_LM >> 1 (command-line tunable, default 25): Allen-Cahn relaxation
 *   is non-conserved (rate ~ L*const at k=0), while Cahn-Hilliard
 *   relaxation is conserved/diffusive (rate ~ M*k^2*curvature -> 0 as
 *   k->0), so c is intrinsically the slow field and eta the fast one --
 *   making L >> M pins eta adiabatically to whichever CALPHAD well c is
 *   closest to at every instant. Naively cranking L up would blow up an
 *   explicit-Euler treatment of the double-well reaction term though
 *   (its curvature w*g''(eta) is O(w) and enters with rate ~L*w), so we
 *   use eigenvalue/convexity splitting: g'(eta) = A0*eta + [g'(eta)-A0*eta]
 *   with A0 = max_eta g''(eta) = 2, and treat the A0*eta part IMPLICITLY
 *   (folded into the spectral denominator alongside kappa_eta*k^2); the
 *   remainder g'(eta)-A0*eta is provably non-stiff (its own derivative is
 *   <=0 everywhere on [0,1]), so it can stay explicit at any L.
 *
 *  MODE "ch": plain (non-KKS) Cahn-Hilliard on a single conserved field.
 *   No eta, no partitioning -- the CALPHAD G(x,T) is simply non-convex
 *   inside the miscibility gap, and that non-convexity alone drives
 *   spinodal decomposition:
 *     dc/dt = div( M grad(mu) ),   mu = dG/dx(c,T) - kappa_c*lap(c)
 *   This is the more standard/simpler way to model a single-phase
 *   miscibility gap; KKS earns its extra complexity only when a second,
 *   structurally distinct phase must be tracked.
 *
 * NUMERICS: semi-implicit Fourier-spectral time stepping (Chen & Shen,
 *   Comp. Phys. Comm. 108 (1998) 147; Zhu, Chen, Shen, Comp. Mater. Sci.
 *   20 (2001)):
 *     chat^{n+1}   = [chat^n   - dt*M*k^2*muhat_bulk^n]
 *                    / (1 + dt*M*kappa_c*k^4)
 *     etahat^{n+1} = [etahat^n - dt*L*(driv_hat^n - w*A0*etahat^n)]
 *                    / (1 + dt*L*(kappa_eta*k^2 + w*A0))          [kks mode]
 *   The stiff gradient-energy (and, for eta, double-well curvature) terms
 *   are treated implicitly (unconditionally stable); the nonlinear bulk
 *   CALPHAD driving forces are treated explicitly.
 *
 * BUILD:
 *   gcc -O3 -march=native -o kks_fecr_spinodal kks_fecr_spinodal.c \
 *       -lfftw3 -lm
 *   (Debian/Ubuntu: sudo apt-get install libfftw3-dev)
 *
 * RUN:
 *   ./kks_fecr_spinodal [mode: kks|ch] [Nx] [Ny] [T_kelvin] [c0] \
 *                       [nsteps] [out_every] [L_over_M (kks only)]
 *   e.g.  ./kks_fecr_spinodal kks 128 128 700 0.45 20000 500 25
 *         ./kks_fecr_spinodal ch  128 128 700 0.45 20000 500
 *
 * OUTPUT: legacy-VTK structured-points files kks_out/c_XXXXXX.vtk (and,
 *   in kks mode, eta_XXXXXX.vtk) that can be opened directly in
 *   ParaView/VisIt, plus a bulk-free-energy-vs-time log
 *   (kks_out/energy.dat) for a sanity check (should be monotonically
 *   non-increasing; it excludes the interfacial gradient-energy term).
 * ===================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <fftw3.h>

/* =====================================================================
 * PART 1: minimal TDB parser  (FUNCTION / PAR statements only)
 * ===================================================================== */
#define MAX_TERMS 32
#define MAX_SEG   4
#define NAME_LEN  32

typedef struct {
    double coeff;
    int    tpow;
    int    has_ln;
    char   ref[NAME_LEN];
} Term;

typedef struct {
    double lo, hi;
    Term terms[MAX_TERMS];
    int nterms;
} Segment;

typedef struct {
    char name[NAME_LEN];
    Segment seg[MAX_SEG];
    int nseg;
} TFunc;

#define MAX_FUNCS 64
static TFunc g_funcs[MAX_FUNCS];
static int g_nfuncs = 0;

typedef struct {
    char type[8];
    char phase[24];
    char key[64];
    int  order;
    char exprtext[256];
} ParEntry;
#define MAX_PAR 64
static ParEntry g_pars[MAX_PAR];
static int g_npar = 0;

static TFunc* find_func(const char* name) {
    for (int i = 0; i < g_nfuncs; i++)
        if (strcmp(g_funcs[i].name, name) == 0) return &g_funcs[i];
    return NULL;
}
static double eval_func(const char* name, double T);

static double eval_terms(Segment* s, double T) {
    double v = 0.0;
    for (int i = 0; i < s->nterms; i++) {
        Term* t = &s->terms[i];
        double f = t->coeff;
        if (t->ref[0]) f *= eval_func(t->ref, T);
        if (t->tpow != 0) f *= pow(T, (double)t->tpow);
        if (t->has_ln) f *= log(T);
        v += f;
    }
    return v;
}
static double eval_func(const char* name, double T) {
    TFunc* fn = find_func(name);
    if (!fn) { fprintf(stderr, "TDB: unknown function %s\n", name); exit(1); }
    for (int i = 0; i < fn->nseg; i++)
        if (T <= fn->seg[i].hi + 1e-9 || i == fn->nseg - 1)
            return eval_terms(&fn->seg[i], T);
    return eval_terms(&fn->seg[fn->nseg-1], T);
}

static void trim(char* s) {
    int n = (int)strlen(s);
    int a = 0; while (a < n && isspace((unsigned char)s[a])) a++;
    int b = n - 1; while (b >= a && isspace((unsigned char)s[b])) b--;
    int len = b - a + 1; if (len < 0) len = 0;
    memmove(s, s + a, len); s[len] = 0;
}
static void strip_ws(char* s) {
    char* w = s;
    for (char* r = s; *r; r++) if (!isspace((unsigned char)*r)) *w++ = *r;
    *w = 0;
}
static int split_terms(const char* expr, char terms[][64], int max_terms) {
    int n = (int)strlen(expr);
    int depth = 0, nterm = 0, start = 0;
    for (int i = 0; i < n; i++) {
        char c = expr[i];
        if (c == '(') depth++;
        if (c == ')') depth--;
        if ((c == '+' || c == '-') && depth == 0 && i > start &&
            expr[i-1] != 'E' && expr[i-1] != 'e') {
            int len = i - start; if (len > 63) len = 63;
            strncpy(terms[nterm], expr + start, len); terms[nterm][len] = 0;
            nterm++; if (nterm >= max_terms) return nterm;
            start = i;
        }
    }
    int len = n - start; if (len > 63) len = 63;
    strncpy(terms[nterm], expr + start, len); terms[nterm][len] = 0;
    nterm++;
    return nterm;
}
static int is_number(const char* s, double* out) {
    char* end; double v = strtod(s, &end);
    if (end != s && *end == 0) { *out = v; return 1; }
    return 0;
}
static void parse_term(const char* raw, Term* t) {
    char buf[64]; strncpy(buf, raw, 63); buf[63] = 0;
    double sign = 1.0; int off = 0;
    if (buf[0] == '+') off = 1;
    else if (buf[0] == '-') { sign = -1.0; off = 1; }
    char body[64]; strcpy(body, buf + off);
    char rep[80]; int rj = 0;
    for (int i = 0; body[i]; i++) {
        if (body[i] == '*' && body[i+1] == '*') { rep[rj++] = '^'; i++; }
        else rep[rj++] = body[i];
    }
    rep[rj] = 0;
    char* saveptr; char factors[8][64]; int nf = 0;
    char* tok = strtok_r(rep, "*", &saveptr);
    while (tok && nf < 8) { strcpy(factors[nf++], tok); tok = strtok_r(NULL, "*", &saveptr); }
    double coeff = 1.0; int tpow = 0; int has_ln = 0; char ref[NAME_LEN]; ref[0] = 0;
    for (int i = 0; i < nf; i++) {
        char* f = factors[i]; double num;
        if (strcmp(f, "T") == 0) tpow += 1;
        else if (f[0] == 'T' && f[1] == '^') {
            char expo[32]; strcpy(expo, f + 2);
            int l = (int)strlen(expo);
            if (l > 0 && expo[0] == '(' && expo[l-1] == ')') { expo[l-1] = 0; memmove(expo, expo+1, l-1); }
            tpow += atoi(expo);
        } else if (strcmp(f, "LN(T)") == 0) has_ln = 1;
        else if (is_number(f, &num)) coeff *= num;
        else strncpy(ref, f, NAME_LEN-1);
    }
    t->coeff = sign * coeff; t->tpow = tpow; t->has_ln = has_ln;
    memset(t->ref,0,NAME_LEN); strncpy(t->ref, ref, NAME_LEN-1);
}
static void parse_expr_into_segment(const char* expr, Segment* seg) {
    char terms[MAX_TERMS][64];
    char clean[512]; strncpy(clean, expr, 511); clean[511]=0; strip_ws(clean);
    int nt = split_terms(clean, terms, MAX_TERMS);
    seg->nterms = 0;
    for (int i = 0; i < nt; i++) {
        if (strlen(terms[i]) == 0) continue;
        parse_term(terms[i], &seg->terms[seg->nterms]); seg->nterms++;
    }
}
static char* read_whole_file(const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open TDB file '%s'\n", path); exit(1); }
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    char* buf = malloc(sz + 1);
    size_t rd = fread(buf, 1, sz, f); buf[rd] = 0; fclose(f);
    return buf;
}
static void parse_function_stmt(const char* stmt) {
    const char* p = stmt;
    while (*p && !isspace((unsigned char)*p)) p++;
    while (isspace((unsigned char)*p)) p++;
    char name[NAME_LEN]; int ni = 0;
    while (*p && !isspace((unsigned char)*p) && ni < NAME_LEN-1) name[ni++] = *p++;
    name[ni] = 0;
    while (isspace((unsigned char)*p)) p++;

    TFunc* fn = &g_funcs[g_nfuncs++];
    memset(fn,0,sizeof(*fn));
    strncpy(fn->name, name, NAME_LEN-1);
    fn->nseg = 0;

    char* endp;
    double low = strtod(p, &endp); p = endp;
    while (*p) {
        while (isspace((unsigned char)*p)) p++;
        const char* estart = p;
        while (*p && *p != ';') p++;
        char expr[512]; int len = (int)(p - estart); if (len > 511) len = 511;
        strncpy(expr, estart, len); expr[len] = 0;
        if (*p == ';') p++;
        while (isspace((unsigned char)*p)) p++;
        double high = strtod(p, &endp); p = endp;
        while (isspace((unsigned char)*p)) p++;
        char flag = *p; if (*p) p++;

        Segment* seg = &fn->seg[fn->nseg++];
        seg->lo = low; seg->hi = high;
        parse_expr_into_segment(expr, seg);

        low = high;
        if (flag != 'Y') break;
        while (isspace((unsigned char)*p)) p++;
    }
}
static void parse_par_stmt(const char* stmt) {
    const char* p = stmt;
    while (*p && !isspace((unsigned char)*p)) p++;
    while (isspace((unsigned char)*p)) p++;
    char type[8]; int ti = 0;
    while (*p && *p != '(' && ti < 7) type[ti++] = *p++;
    type[ti] = 0; trim(type);
    if (*p != '(') return;
    p++;
    const char* kstart = p;
    while (*p && *p != ')') p++;
    char key[64]; int klen = (int)(p - kstart); if (klen > 63) klen = 63;
    strncpy(key, kstart, klen); key[klen] = 0;
    if (*p == ')') p++;

    char phase[24] = "", constit[40] = ""; int order = 0;
    char keycopy[64]; strncpy(keycopy,key,63); keycopy[63]=0;
    char* semi = strchr(keycopy, ';');
    if (semi) { order = atoi(semi + 1); *semi = 0; }
    char* comma = strchr(keycopy, ',');
    if (comma) {
        int plen = (int)(comma - keycopy);
        strncpy(phase, keycopy, plen); phase[plen] = 0;
        strcpy(constit, comma + 1);
    } else strcpy(phase, keycopy);

    while (*p == ',' || isspace((unsigned char)*p)) p++;
    const char* estart = p;
    while (*p && *p != ';') p++;
    char expr[256]; int elen = (int)(p - estart); if (elen > 255) elen = 255;
    strncpy(expr, estart, elen); expr[elen] = 0;

    ParEntry* pe = &g_pars[g_npar++];
    memset(pe,0,sizeof(*pe));
    strcpy(pe->type, type); strcpy(pe->phase, phase);
    strcpy(pe->key, constit); pe->order = order;
    strcpy(pe->exprtext, expr);
}
static void parse_tdb(const char* path) {
    char* raw = read_whole_file(path);
    char* joined = malloc(strlen(raw) * 2 + 16);
    int jl = 0;
    char* line = strtok(raw, "\n");
    while (line) {
        char tmp[1024]; strncpy(tmp, line, 1023); tmp[1023] = 0; trim(tmp);
        if (tmp[0] != '$' && tmp[0] != 0) {
            int l = (int)strlen(tmp);
            memcpy(joined + jl, tmp, l); jl += l;
            joined[jl++] = ' ';
        }
        line = strtok(NULL, "\n");
    }
    joined[jl] = 0;
    char* stmt = strtok(joined, "!");
    while (stmt) {
        char s[2048]; strncpy(s, stmt, 2047); s[2047]=0; trim(s);
        if (strlen(s) > 0) {
            if (strncasecmp(s, "FUNCTION", 8) == 0) parse_function_stmt(s);
            else if (strncasecmp(s, "PAR", 3) == 0) parse_par_stmt(s);
        }
        stmt = strtok(NULL, "!");
    }
    free(raw); free(joined);
}
static ParEntry* find_par(const char* type, const char* key, int order) {
    for (int i = 0; i < g_npar; i++)
        if (strcmp(g_pars[i].type, type) == 0 &&
            strcmp(g_pars[i].key, key) == 0 &&
            g_pars[i].order == order) return &g_pars[i];
    return NULL;
}
static double eval_par_expr(ParEntry* pe, double T) {
    if (!pe) return 0.0;
    Segment seg; parse_expr_into_segment(pe->exprtext, &seg);
    return eval_terms(&seg, T);
}

/* =====================================================================
 * PART 2: CALPHAD molar Gibbs energy of BCC_A2 Fe-Cr (from parsed TDB)
 * ===================================================================== */
#define RGAS 8.31451

static double magnetic_g(double tau, double p) {
    double D = 518.0/1125.0 + (11692.0/15975.0) * (1.0/p - 1.0);
    if (tau <= 1.0) {
        double term = 79.0*pow(tau, -1.0)/(140.0*p)
                    + (474.0/497.0)*(1.0/p - 1.0) *
                      (pow(tau,3)/6.0 + pow(tau,9)/135.0 + pow(tau,15)/600.0);
        return 1.0 - term / D;
    } else {
        double term = pow(tau,-5.0)/10.0 + pow(tau,-15.0)/315.0 + pow(tau,-25.0)/1500.0;
        return -term / D;
    }
}

/* cached PAR pointers, resolved once after parsing */
static ParEntry *pL0, *pTcFe, *pTcCr, *pBmFe, *pBmCr, *pTc0, *pTc1, *pBm0;

static void calphad_init(const char* tdb_path) {
    parse_tdb(tdb_path);
    pL0   = find_par("L", "CR,FE:VA", 0);
    pTcFe = find_par("TC","FE:VA",0);
    pTcCr = find_par("TC","CR:VA",0);
    pBmFe = find_par("BM","FE:VA",0);
    pBmCr = find_par("BM","CR:VA",0);
    pTc0  = find_par("TC","CR,FE:VA",0);
    pTc1  = find_par("TC","CR,FE:VA",1);
    pBm0  = find_par("BM","CR,FE:VA",0);
    if (!pL0 || !pTcFe || !pTcCr || !pBmFe || !pBmCr || !pTc0 || !pTc1 || !pBm0) {
        fprintf(stderr, "TDB: required Fe-Cr BCC_A2 parameters missing from %s\n", tdb_path);
        exit(1);
    }
}

/* molar Gibbs energy [J/mol], x = mole fraction Cr, clamped away from 0,1 */
static double G_bcc(double x, double T) {
    if (x < 1e-6) x = 1e-6;
    if (x > 1 - 1e-6) x = 1 - 1e-6;
    double xFe = 1.0 - x, xCr = x;

    double g_ref = xFe*eval_func("GHSERFE", T) + xCr*eval_func("GHSERCR", T);
    double g_ideal = RGAS*T*(xCr*log(xCr) + xFe*log(xFe));
    double L0 = eval_par_expr(pL0, T);
    double g_excess = xCr*xFe*L0;

    double tcFe = eval_par_expr(pTcFe, T), tcCr = eval_par_expr(pTcCr, T);
    double bmFe = eval_par_expr(pBmFe, T), bmCr = eval_par_expr(pBmCr, T);
    double tc0 = eval_par_expr(pTc0, T), tc1 = eval_par_expr(pTc1, T);
    double bm0 = eval_par_expr(pBm0, T);

    double Tc = xFe*tcFe + xCr*tcCr + xCr*xFe*(tc0 + tc1*(xCr - xFe));
    double beta = xFe*bmFe + xCr*bmCr + xCr*xFe*bm0;
    if (Tc < 1.0) Tc = 1.0;
    if (beta < 1e-6) beta = 1e-6;
    double g_mag = RGAS*T*log(beta + 1.0) * magnetic_g(T/Tc, 0.4);

    return g_ref + g_ideal + g_excess + g_mag;
}

/* first/second derivative of G wrt x at fixed T (central finite differences) */
static double dG_dx(double x, double T) {
    double h = 1e-5;
    if (x < h) x = h; if (x > 1-h) x = 1-h;
    return (G_bcc(x+h,T) - G_bcc(x-h,T)) / (2*h);
}
static double d2G_dx2(double x, double T) {
    double h = 2e-4;
    if (x < h) x = h; if (x > 1-h) x = 1-h;
    return (G_bcc(x+h,T) - 2*G_bcc(x,T) + G_bcc(x-h,T)) / (h*h);
}

/* =====================================================================
 * PART 3: KKS local equilibrium partition solver (2x2 Newton)
 * =====================================================================
 * Given overall composition c and phase indicator eta, find (cA,cB) with
 *   c = h(eta) cB + (1-h(eta)) cA
 *   dG/dx(cA) = dG/dx(cB)
 * Both "phases" A and B use the SAME CALPHAD G(x,T) (single physical BCC
 * phase); only their local compositions differ -- this is exactly how KKS
 * represents a single-phase miscibility gap / spinodal decomposition.
 */
static double h_interp(double e) { return e*e*e*(6.0*e*e - 15.0*e + 10.0); }
static double hp_interp(double e) { return 30.0*e*e*(1.0-e)*(1.0-e); }
static double g_doublewell(double e) { return e*e*(1.0-e)*(1.0-e); }
static double gp_doublewell(double e) { return 2.0*e*(1.0-e)*(1.0-2.0*e); }

static void kks_partition(double c, double eta, double T,
                           double* cA_io, double* cB_io) {
    double h = h_interp(eta);
    double cA = *cA_io, cB = *cB_io;
    /* clamp initial guesses */
    if (cA < 0.001) cA = 0.001; if (cA > 0.999) cA = 0.999;
    if (cB < 0.001) cB = 0.001; if (cB > 0.999) cB = 0.999;

    for (int it = 0; it < 25; it++) {
        double R1 = h*cB + (1.0-h)*cA - c;
        double muA = dG_dx(cA, T), muB = dG_dx(cB, T);
        double R2 = muB - muA;
        if (fabs(R1) < 1e-10 && fabs(R2) < 1e-4) break;

        double J11 = (1.0-h), J12 = h;
        double J21 = -d2G_dx2(cA, T), J22 = d2G_dx2(cB, T);
        double det = J11*J22 - J12*J21;
        if (fabs(det) < 1e-12) break;
        double dA = ( J22*(-R1) - J12*(-R2)) / det;
        double dB = (-J21*(-R1) + J11*(-R2)) / det;

        /* damped update, keep inside (0,1) */
        double step = 1.0;
        double nA = cA + step*dA, nB = cB + step*dB;
        while ((nA <= 1e-4 || nA >= 1-1e-4 || nB <= 1e-4 || nB >= 1-1e-4) && step > 1e-4) {
            step *= 0.5; nA = cA + step*dA; nB = cB + step*dB;
        }
        cA = nA; cB = nB;
    }
    *cA_io = cA; *cB_io = cB;
}

/* =====================================================================
 * PART 4: semi-implicit Fourier-spectral solver -- "kks" or "ch" mode
 * ===================================================================== */
typedef struct { int Nx, Ny; double dx, dy; } Grid;
typedef enum { MODE_KKS, MODE_CH } Mode;

static void write_vtk(const char* fname, const double* field, Grid g, const char* label) {
    FILE* f = fopen(fname, "w");
    if (!f) { fprintf(stderr, "cannot open %s\n", fname); return; }
    fprintf(f, "# vtk DataFile Version 3.0\nKKS %s\nASCII\nDATASET STRUCTURED_POINTS\n", label);
    fprintf(f, "DIMENSIONS %d %d 1\n", g.Nx, g.Ny);
    fprintf(f, "ORIGIN 0 0 0\nSPACING %g %g 1\n", g.dx, g.dy);
    fprintf(f, "POINT_DATA %d\nSCALARS %s double 1\nLOOKUP_TABLE default\n", g.Nx*g.Ny, label);
    for (int j = 0; j < g.Ny; j++)
        for (int i = 0; i < g.Nx; i++)
            fprintf(f, "%.6f\n", field[j*g.Nx+i]);
    fclose(f);
}

int main(int argc, char** argv) {
    setvbuf(stdout, NULL, _IOLBF, 0);

    Mode mode = MODE_KKS;
    if (argc > 1) {
        if      (strcasecmp(argv[1], "ch")  == 0) mode = MODE_CH;
        else if (strcasecmp(argv[1], "kks") == 0) mode = MODE_KKS;
        else { fprintf(stderr, "unknown mode '%s' (use kks|ch)\n", argv[1]); return 1; }
    }
    int Nx        = argc>2 ? atoi(argv[2]) : 128;
    int Ny        = argc>3 ? atoi(argv[3]) : 128;
    double T      = argc>4 ? atof(argv[4]) : 700.0;   /* Kelvin, inside spinodal (~<900K) */
    double c0     = argc>5 ? atof(argv[5]) : 0.45;    /* mean Cr mole fraction */
    int nsteps    = argc>6 ? atoi(argv[6]) : 20000;
    int out_every = argc>7 ? atoi(argv[7]) : 500;
    double R_LM   = argc>8 ? atof(argv[8]) : 25.0;    /* L/M ratio, kks mode only */
    const char* tdb_path = "fecr_bcc.tdb";

    calphad_init(tdb_path);
    printf("Loaded CALPHAD Fe-Cr BCC_A2 parameters from %s\n", tdb_path); fflush(stdout);
    printf("mode=%s  T=%.1fK  c0(Cr)=%.3f  grid=%dx%d  nsteps=%d%s\n",
           mode==MODE_KKS ? "kks" : "ch", T, c0, Nx, Ny, nsteps,
           mode==MODE_KKS ? "" : "");
    if (mode == MODE_KKS) printf("L/M = %.1f  (eta relaxes ~%.0fx faster than c)\n", R_LM, R_LM);

    double d2g0 = d2G_dx2(c0, T);
    printf("d2G/dx2(c0,T) = %.1f J/mol  (%s => spinodally unstable start point)\n",
           d2g0, d2g0 < 0 ? "negative" : "positive");

    mkdir("kks_out", 0755);

    Grid grid = { Nx, Ny, 1.0, 1.0 };
    int N = Nx*Ny;

    /* ---- nondimensionalization: energies in units of E0 = R*T ---- */
    double E0 = RGAS*T;

    /* ---- phase-field model parameters (demo values; calibrate against
     *      real interfacial energy/mobility data for quantitative work) */
    double M       = 1.0;          /* atomic mobility (nondimensional)          */
    double L       = R_LM * M;     /* eta relaxation coeff: L >> M (kks only)   */
    double kappa_c = 2.0;          /* gradient energy coeff for c (grid units)  */
    double kappa_e = 2.0;          /* gradient energy coeff for eta             */
    double w       = 2.0;          /* double-well barrier height (units of E0)  */
    double dt      = 0.5;
    const double A0 = 2.0;         /* convexity-splitting constant: max g''(eta) on [0,1] */

    /* ---- fields (eta/cA/cB only used in kks mode) ---- */
    double *c   = fftw_malloc(sizeof(double)*N);
    double *eta = NULL, *cA = NULL, *cB = NULL, *eta_driv = NULL;
    double *mu_bulk  = fftw_malloc(sizeof(double)*N);   /* dG/dx(c or cA)/E0 */

    if (mode == MODE_KKS) {
        eta = fftw_malloc(sizeof(double)*N);
        cA  = fftw_malloc(sizeof(double)*N);
        cB  = fftw_malloc(sizeof(double)*N);
        eta_driv = fftw_malloc(sizeof(double)*N);
    }

    srand(12345);
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            int idx = j*Nx+i;
            double r1 = ((double)rand()/RAND_MAX - 0.5);
            c[idx] = c0 + 0.01*r1;
            if (mode == MODE_KKS) {
                /* Seed eta from c's OWN fluctuation (not independent noise): with
                 * L >> M, eta saturates to 0/1 almost immediately -- essentially
                 * before c has developed any real structure of its own (see
                 * discussion in header). If eta started from independent random
                 * noise it would lock onto ITS noise pattern, uncorrelated with
                 * c. Keying the seed to sign(c-c0) instead gives it something
                 * physically meaningful to lock onto from step zero. */
                eta[idx] = 0.5 + 0.49*tanh((c[idx] - c0) * 300.0);
                cA[idx] = c0 - 0.05; cB[idx] = c0 + 0.05;
            }
        }

    /* ---- FFTW plans (real-to-complex 2D) ---- */
    int Nyh = Ny/2 + 1;
    fftw_complex *chat  = fftw_malloc(sizeof(fftw_complex)*Nx*Nyh);
    fftw_complex *muhat = fftw_malloc(sizeof(fftw_complex)*Nx*Nyh);
    fftw_plan p_c_fwd  = fftw_plan_dft_r2c_2d(Nx, Ny, c,       chat, FFTW_ESTIMATE);
    fftw_plan p_mu_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, mu_bulk, muhat,FFTW_ESTIMATE);
    fftw_plan p_c_bwd  = fftw_plan_dft_c2r_2d(Nx, Ny, chat, c, FFTW_ESTIMATE);

    fftw_complex *etahat = NULL, *drvhat = NULL;
    fftw_plan p_eta_fwd = 0, p_drv_fwd = 0, p_eta_bwd = 0;
    if (mode == MODE_KKS) {
        etahat = fftw_malloc(sizeof(fftw_complex)*Nx*Nyh);
        drvhat = fftw_malloc(sizeof(fftw_complex)*Nx*Nyh);
        p_eta_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, eta,      etahat, FFTW_ESTIMATE);
        p_drv_fwd = fftw_plan_dft_r2c_2d(Nx, Ny, eta_driv, drvhat, FFTW_ESTIMATE);
        p_eta_bwd = fftw_plan_dft_c2r_2d(Nx, Ny, etahat, eta, FFTW_ESTIMATE);
    }

    /* wavenumbers */
    double *kx2 = malloc(sizeof(double)*Nx);
    double *ky2 = malloc(sizeof(double)*Nyh);
    for (int i = 0; i < Nx; i++) {
        int ii = (i <= Nx/2) ? i : i - Nx;
        double k = 2.0*M_PI*ii/(Nx*grid.dx);
        kx2[i] = k*k;
    }
    for (int j = 0; j < Nyh; j++) {
        double k = 2.0*M_PI*j/(Ny*grid.dy);
        ky2[j] = k*k;
    }

    FILE* efile = fopen("kks_out/energy.dat", "w");
    fprintf(efile, "# step  time  bulk_free_energy[J/mol]  mean_c\n");

    clock_t t_start = clock();
    for (int step = 0; step <= nsteps; step++) {

        double F_total = 0.0, mean_c = 0.0;

        if (mode == MODE_KKS) {
            /* 1a. local KKS equilibrium partition at every grid point */
            for (int idx = 0; idx < N; idx++) {
                kks_partition(c[idx], eta[idx], T, &cA[idx], &cB[idx]);
                double hh = h_interp(eta[idx]);
                double mu = dG_dx(cA[idx], T);
                mu_bulk[idx] = mu / E0;
                double GA = G_bcc(cA[idx], T), GB = G_bcc(cB[idx], T);
                double gp = gp_doublewell(eta[idx]);
                /* convexity split: explicit remainder = w*(g'(eta) - A0*eta) + coupling term */
                eta_driv[idx] = ( hp_interp(eta[idx])*(GB - GA - (cB[idx]-cA[idx])*mu)
                                   + w*E0*(gp - A0*eta[idx]) ) / E0;
                if (step % out_every == 0) {
                    F_total += (hh*GB + (1.0-hh)*GA + w*E0*g_doublewell(eta[idx]));
                    mean_c += c[idx];
                }
            }
        } else {
            /* 1b. plain Cahn-Hilliard: mu = dG/dx(c,T) directly, no partition */
            for (int idx = 0; idx < N; idx++) {
                mu_bulk[idx] = dG_dx(c[idx], T) / E0;
                if (step % out_every == 0) {
                    F_total += G_bcc(c[idx], T);
                    mean_c  += c[idx];
                }
            }
        }

        if (step % out_every == 0) {
            F_total /= N; mean_c /= N;
            printf("step %6d  t=%8.3f  <F_bulk>=%12.4f J/mol  <c>=%.5f\n",
                   step, step*dt, F_total, mean_c);
            fprintf(efile, "%d %g %g %g\n", step, step*dt, F_total, mean_c);
            fflush(efile);

            char fn1[128]; snprintf(fn1, sizeof(fn1), "kks_out/c_%06d.vtk", step);
            write_vtk(fn1, c, grid, "Cr_fraction");
            if (mode == MODE_KKS) {
                char fn2[128]; snprintf(fn2, sizeof(fn2), "kks_out/eta_%06d.vtk", step);
                write_vtk(fn2, eta, grid, "eta");
            }
        }
        if (step == nsteps) break;

        /* 2. forward transforms */
        fftw_execute(p_c_fwd);
        fftw_execute(p_mu_fwd);
        if (mode == MODE_KKS) { fftw_execute(p_eta_fwd); fftw_execute(p_drv_fwd); }

        /* 3. semi-implicit spectral update */
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Nyh; j++) {
                int idx = i*Nyh + j;
                double k2 = kx2[i] + ky2[j];
                double k4 = k2*k2;

                double denom_c = 1.0 + dt*M*kappa_c*k4;
                chat[idx][0] = (chat[idx][0] - dt*M*k2*muhat[idx][0]) / denom_c;
                chat[idx][1] = (chat[idx][1] - dt*M*k2*muhat[idx][1]) / denom_c;

                if (mode == MODE_KKS) {
                    /* convexity-split implicit denominator: gradient term + w*A0 */
                    double denom_e = 1.0 + dt*L*(kappa_e*k2 + w*A0);
                    etahat[idx][0] = (etahat[idx][0] - dt*L*drvhat[idx][0]) / denom_e;
                    etahat[idx][1] = (etahat[idx][1] - dt*L*drvhat[idx][1]) / denom_e;
                }
            }
        }

        /* 4. inverse transforms (+FFTW normalization) */
        fftw_execute(p_c_bwd);
        if (mode == MODE_KKS) fftw_execute(p_eta_bwd);
        double norm = 1.0/(Nx*Ny);
        for (int idx = 0; idx < N; idx++) {
            c[idx] *= norm;
            if (c[idx] < 0.002) c[idx] = 0.002;
            if (c[idx] > 0.998) c[idx] = 0.998;
            if (mode == MODE_KKS) {
                eta[idx] *= norm;
                if (eta[idx] < -0.05) eta[idx] = -0.05;
                if (eta[idx] >  1.05) eta[idx] =  1.05;
            }
        }
    }
    double elapsed = (double)(clock()-t_start)/CLOCKS_PER_SEC;
    printf("Done. %d steps in %.2f s (%.3f ms/step). VTK output in kks_out/\n",
           nsteps, elapsed, 1000.0*elapsed/nsteps);
    fclose(efile);

    fftw_destroy_plan(p_c_fwd); fftw_destroy_plan(p_mu_fwd); fftw_destroy_plan(p_c_bwd);
    fftw_free(chat); fftw_free(muhat);
    fftw_free(c); fftw_free(mu_bulk);
    if (mode == MODE_KKS) {
        fftw_destroy_plan(p_eta_fwd); fftw_destroy_plan(p_drv_fwd); fftw_destroy_plan(p_eta_bwd);
        fftw_free(etahat); fftw_free(drvhat);
        fftw_free(eta); fftw_free(cA); fftw_free(cB); fftw_free(eta_driv);
    }
    free(kx2); free(ky2);
    return 0;
}
