// ============================================================================
//  poly_calc.cpp  --  一元多项式的带余除法 与 最大公因数 (GCD)
//
//  说明：
//    * 系数在有理数域 Q 上运算，支持 整数 / 分数 / 小数 形式输入。
//    * 带余除法： A(x) = B(x) * Q(x) + R(x)，且 deg R < deg B。
//    * 最大公因数：用辗转相除法（欧几里得算法）求 gcd，
//      结果归一化为「首一多项式」(monic)，并同时给出
//      「整数本原形式」（消去分母、提出内容、首项为正）。
//    * 注意：域 Q 上 gcd 只在「相差一个非零常数倍」的意义下唯一，
//      所以程序统一输出首一形式作为代表。
//
//  编译（MinGW / GCC）：
//      g++ -std=c++11 -O2 -o poly_calc poly_calc.cpp
//  编译（MSVC）：
//      cl /std:c++14 /utf-8 /EHsc poly_calc.cpp
//  运行：
//      poly_calc
// ============================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cctype>
#include <climits>

#ifdef _WIN32
#include <windows.h>
#endif

// ---------------------------------------------------------------------------
//  宽整数：中间乘加尽量用 128 位，降低溢出风险
// ---------------------------------------------------------------------------
#if defined(__SIZEOF_INT128__)
typedef __int128 wide_t;
#define HAS_WIDE 1
#else
typedef long long wide_t;
#define HAS_WIDE 0
#endif

static bool g_overflow = false;   // 一旦发生 long long 溢出就置位

static wide_t wabs(wide_t x) { return x < 0 ? -x : x; }

static wide_t wgcd(wide_t a, wide_t b) {
    a = wabs(a); b = wabs(b);
    while (b != 0) { wide_t t = a % b; a = b; b = t; }
    return a;
}

static long long toLL(wide_t x) {
#if HAS_WIDE
    if (x > (wide_t)LLONG_MAX || x < (wide_t)LLONG_MIN) {
        g_overflow = true;
        return (long long)(x < 0 ? LLONG_MIN : LLONG_MAX);
    }
#endif
    return (long long)x;
}

// ---------------------------------------------------------------------------
//  有理数（分数）：p/q，恒约分，分母恒为正
// ---------------------------------------------------------------------------
struct Frac {
    long long p;   // 分子
    long long q;   // 分母 (> 0)

    Frac() : p(0), q(1) {}
    explicit Frac(long long v) : p(0), q(1) { set((wide_t)v, (wide_t)1); }
    Frac(long long pp, long long qq) : p(0), q(1) { set((wide_t)pp, (wide_t)qq); }

    void set(wide_t pp, wide_t qq) {
        if (qq == 0) { pp = 0; qq = 1; }         // 保护：不应发生
        if (qq < 0) { pp = -pp; qq = -qq; }
        if (pp == 0) { p = 0; q = 1; return; }
        wide_t g = wgcd(pp, qq);
        if (g > 1) { pp /= g; qq /= g; }
        p = toLL(pp); q = toLL(qq);
    }

    bool isZero() const { return p == 0; }
};

static Frac addF(const Frac& a, const Frac& b) {
    Frac r; r.set((wide_t)a.p * b.q + (wide_t)b.p * a.q, (wide_t)a.q * b.q); return r;
}
static Frac subF(const Frac& a, const Frac& b) {
    Frac r; r.set((wide_t)a.p * b.q - (wide_t)b.p * a.q, (wide_t)a.q * b.q); return r;
}
static Frac mulF(const Frac& a, const Frac& b) {
    // 交叉约分：先约掉公因子，再相乘，减小中间结果
    wide_t n1 = a.p, d1 = a.q, n2 = b.p, d2 = b.q;
    wide_t g1 = wgcd(n1, d2);
    if (g1 > 1) { n1 /= g1; d2 /= g1; }
    wide_t g2 = wgcd(n2, d1);
    if (g2 > 1) { n2 /= g2; d1 /= g2; }
    Frac r; r.set(n1 * n2, d1 * d2); return r;
}
static Frac divF(const Frac& a, const Frac& b) {
    // a / b = (a.p * b.q) / (a.q * b.p)
    wide_t num = a.p, den = a.q, bn = b.p, bd = b.q;
    wide_t g1 = wgcd(num, bn);
    if (g1 > 1) { num /= g1; bn /= g1; }
    wide_t g2 = wgcd(bd, den);
    if (g2 > 1) { bd /= g2; den /= g2; }
    Frac r; r.set(num * bd, den * bn); return r;
}

// ---------------------------------------------------------------------------
//  多项式：c[i] 是 x^i 的系数；空 vector 表示零多项式
// ---------------------------------------------------------------------------
typedef std::vector<Frac> Poly;

static void trim(Poly& a) {
    while (!a.empty() && a.back().isZero()) a.pop_back();
}
static int degOf(const Poly& a) { return (int)a.size() - 1; }   // 零多项式为 -1

static Poly addPoly(const Poly& a, const Poly& b) {
    Poly c(std::max(a.size(), b.size()), Frac(0));
    for (size_t i = 0; i < a.size(); ++i) c[i] = addF(c[i], a[i]);
    for (size_t i = 0; i < b.size(); ++i) c[i] = addF(c[i], b[i]);
    trim(c); return c;
}
static Poly mulPoly(const Poly& a, const Poly& b) {
    if (a.empty() || b.empty()) return Poly();
    Poly c(a.size() + b.size() - 1, Frac(0));
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j)
            c[i + j] = addF(c[i + j], mulF(a[i], b[j]));
    trim(c); return c;
}
static Poly subPoly(const Poly& a, const Poly& b) {
    Poly c(std::max(a.size(), b.size()), Frac(0));
    for (size_t i = 0; i < a.size(); ++i) c[i] = addF(c[i], a[i]);
    for (size_t i = 0; i < b.size(); ++i) c[i] = subF(c[i], b[i]);
    trim(c); return c;
}
static Poly scalePoly(const Poly& a, const Frac& k) {
    if (k.isZero()) return Poly();
    Poly c = a;
    for (size_t i = 0; i < c.size(); ++i) c[i] = mulF(c[i], k);
    trim(c); return c;
}
static bool eqPoly(const Poly& a, const Poly& b) {
    Poly x = a, y = b; trim(x); trim(y);
    if (x.size() != y.size()) return false;
    for (size_t i = 0; i < x.size(); ++i)
        if (!(x[i].p == y[i].p && x[i].q == y[i].q)) return false;
    return true;
}

// 带余除法： a = b * q + r ,  deg r < deg b 。 返回 false 表示 b 为零多项式。
static bool polyDivMod(const Poly& a, const Poly& b, Poly& q, Poly& r) {
    if (b.empty()) return false;
    Poly rem = a;
    int db = (int)b.size() - 1;
    int da = (int)a.size() - 1;
    q.assign((size_t)std::max(0, da - db + 1), Frac(0));

    for (int i = da; i >= db; --i) {
        if (i >= (int)rem.size() || rem[i].isZero()) continue;
        Frac factor = divF(rem[i], b[db]);        // 消去当前最高次项
        int shift = i - db;
        q[shift] = factor;
        for (int j = 0; j <= db; ++j)
            rem[shift + j] = subF(rem[shift + j], mulF(factor, b[j]));
    }
    trim(q);
    trim(rem);
    r = rem;
    return true;
}

// 除以首项系数 → 首一多项式
static Poly monic(const Poly& a) {
    Poly r = a; trim(r);
    if (r.empty()) return r;
    Frac lead = r.back();
    for (size_t i = 0; i < r.size(); ++i) r[i] = divF(r[i], lead);
    return r;
}

// 整数本原形式：乘上分母的最小公倍数，再除以分子们的最大公约数，首项取正
static Poly primitiveInteger(const Poly& a) {
    Poly r = a; trim(r);
    if (r.empty()) return r;

    wide_t l = 1;                                  // 分母最小公倍数
    for (size_t i = 0; i < r.size(); ++i) {
        wide_t q = r[i].q;
        wide_t g = wgcd(l, q);
        l = l / (g ? g : 1) * q;
    }
    std::vector<wide_t> num(r.size());
    for (size_t i = 0; i < r.size(); ++i)
        num[i] = (wide_t)r[i].p * (l / r[i].q);

    wide_t g = 0;                                  // 分子最大公约数
    for (size_t i = 0; i < num.size(); ++i) g = wgcd(g, num[i]);
    if (g == 0) g = 1;
    if (num.back() < 0) g = -g;                    // 让首项为正

    Poly out;
    for (size_t i = 0; i < num.size(); ++i) out.push_back(Frac(toLL(num[i] / g)));
    trim(out);
    return out;
}

// 前置声明（供 polyGcd 打印中间过程时使用）
static std::string polyToStr(const Poly& a);

// 辗转相除法求 gcd，结果取首一形式
static Poly polyGcd(Poly a, Poly b, std::vector<std::string>* log) {
    trim(a); trim(b);
    int step = 0;
    while (!b.empty()) {
        Poly q, r;
        polyDivMod(a, b, q, r);
        ++step;
        if (log) {
            std::ostringstream os;
            os << "第 " << step << " 步： A = B * Q + R\n";
            os << "      A = " << polyToStr(a) << "\n";   // 前置声明见下
            os << "      B = " << polyToStr(b) << "\n";
            os << "      Q = " << polyToStr(q) << "\n";
            os << "      R = " << polyToStr(r) << "\n";
            log->push_back(os.str());
        }
        a = b;
        b = r;
    }
    return monic(a);
}
// 扩展欧几里得算法：求 u, v 使 u(x)*f(x) + v(x)*g(x) = d(x) = gcd(f, g)（首一）
// 原理：对余式序列 r[i] = r[i-2] - q[i]*r[i-1] 同步维护
//       r[i] = s[i]*f + t[i]*g，最后把 r[k] 归一化为首一即可。
static void polyExtGcd(const Poly& f, const Poly& g, Poly& u, Poly& v, Poly& d) {
    Poly r0 = f, r1 = g;
    Poly s0(1, Frac(1)), s1;        // s0 = 1,  s1 = 0
    Poly t0, t1(1, Frac(1));        // t0 = 0,  t1 = 1
    trim(r0); trim(r1);

    while (!r1.empty()) {
        Poly q, r;
        polyDivMod(r0, r1, q, r);   // r0 = r1 * q + r
        Poly s = subPoly(s0, mulPoly(q, s1));
        Poly t = subPoly(t0, mulPoly(q, t1));
        r0 = r1; r1 = r;
        s0 = s1; s1 = s;
        t0 = t1; t1 = t;
    }

    u.clear(); v.clear(); d.clear();
    if (r0.empty()) return;                 // f = g = 0 的特例
    Frac inv = divF(Frac(1), r0.back());    // 除以首项系数，使 d 成为首一多项式
    u = scalePoly(s0, inv);
    v = scalePoly(t0, inv);
    d = scalePoly(r0, inv);
}

// ---------------------------------------------------------------------------
//  输出 / 输入
// ---------------------------------------------------------------------------
static std::string llToStr(long long v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    unsigned long long u = neg ? (unsigned long long)(-(v + 1)) + 1ULL : (unsigned long long)v;
    std::string s;
    while (u > 0) { s += (char)('0' + (int)(u % 10)); u /= 10; }
    if (neg) s += '-';
    std::reverse(s.begin(), s.end());
    return s;
}
static std::string fracToStr(const Frac& f) {
    if (f.q == 1) return llToStr(f.p);
    return llToStr(f.p) + "/" + llToStr(f.q);
}

static std::string polyToStr(const Poly& a) {
    if (a.empty()) return "0";
    std::string s;
    for (int i = (int)a.size() - 1; i >= 0; --i) {
        if (a[i].isZero()) continue;
        Frac f = a[i];
        bool neg = f.p < 0;
        Frac mag = neg ? Frac(-f.p, f.q) : f;            // 取绝对值
        if (s.empty()) { if (neg) s += "-"; }
        else s += (neg ? " - " : " + ");

        if (i == 0) {
            s += fracToStr(mag);
        } else {
            bool isOne = (mag.p == 1 && mag.q == 1);
            if (!isOne) {
                if (mag.q == 1) s += llToStr(mag.p);
                else            s += "(" + fracToStr(mag) + ")";
            }
            s += "x";
            if (i > 1) { s += "^"; s += llToStr(i); }
        }
    }
    return s.empty() ? "0" : s;
}

static bool parseIntStr(const std::string& t, long long& out) {
    if (t.empty()) return false;
    size_t k = 0; bool neg = false;
    if (t[k] == '+' || t[k] == '-') { neg = (t[k] == '-'); ++k; }
    if (k >= t.size()) return false;
    unsigned long long v = 0;
    for (; k < t.size(); ++k) {
        if (!isdigit((unsigned char)t[k])) return false;
        v = v * 10ULL + (unsigned long long)(t[k] - '0');
        if (v > 0x7fffffffffffffffULL) return false;
    }
    out = neg ? -(long long)v : (long long)v;
    return true;
}

// 支持 "3" , "-5" , "3/4" , "1.25" , "-0.5"
static bool parseFrac(const std::string& tk, Frac& out) {
    if (tk.empty()) return false;

    size_t slash = tk.find('/');
    if (slash != std::string::npos) {
        long long a, b;
        if (!parseIntStr(tk.substr(0, slash), a)) return false;
        if (!parseIntStr(tk.substr(slash + 1), b)) return false;
        if (b == 0) return false;
        out = Frac(a, b);
        return true;
    }

    size_t dot = tk.find('.');
    if (dot != std::string::npos) {
        std::string digits = tk.substr(0, dot) + tk.substr(dot + 1);
        bool neg = false;
        if (!digits.empty() && (digits[0] == '+' || digits[0] == '-')) {
            neg = (digits[0] == '-');
            digits = digits.substr(1);
        }
        if (digits.empty()) return false;
        long long num = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            if (!isdigit((unsigned char)digits[i])) return false;
            num = num * 10 + (digits[i] - '0');
        }
        long long den = 1;
        for (size_t i = dot + 1; i < tk.size(); ++i) den *= 10;
        if (neg) num = -num;
        out = Frac(num, den);
        return true;
    }

    long long a;
    if (!parseIntStr(tk, a)) return false;
    out = Frac(a);
    return true;
}

// 一行空格分隔的系数，按「从高次到低次」读入
static bool parsePolyLine(const std::string& line, Poly& out, std::string& err) {
    std::istringstream iss(line);
    std::string tk;
    std::vector<Frac> cs;
    while (iss >> tk) {
        Frac f;
        if (!parseFrac(tk, f)) { err = "无法识别的系数: \"" + tk + "\""; return false; }
        cs.push_back(f);
    }
    std::reverse(cs.begin(), cs.end());   // c[i] 对应 x^i
    out = cs;
    trim(out);
    return true;
}

// ---------------------------------------------------------------------------
int main() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    std::cout << "========================================================\n";
    std::cout << "  一元多项式 带余除法 + 最大公因数 计算器 (系数可为分数)\n";
    std::cout << "========================================================\n";
    std::cout << "输入方式：按「从高次到低次」输入全部系数，空格分隔。\n";
    std::cout << "  例:  1 -2 3/4 5   表示  x^3 - 2x^2 + (3/4)x + 5\n";
    std::cout << "  系数支持 整数 / 分数(3/4) / 小数(0.75)。\n";
    std::cout << "  只输入一个 0 或直接回车，表示零多项式。\n\n";

    std::string lineA, lineB;
    for (;;) {
        std::cout << "请输入多项式 A 的系数 > ";
        if (!std::getline(std::cin, lineA)) break;
        Poly A; std::string err;
        if (!parsePolyLine(lineA, A, err)) {
            std::cout << "[错误] A 输入有误：" << err << "\n\n";
            continue;
        }

        std::cout << "请输入多项式 B 的系数 > ";
        if (!std::getline(std::cin, lineB)) break;
        Poly B;
        if (!parsePolyLine(lineB, B, err)) {
            std::cout << "[错误] B 输入有误：" << err << "\n\n";
            continue;
        }

        std::cout << "\n--------------------------------------------------------\n";
        std::cout << "A(x) = " << polyToStr(A)
                  << "      [deg A = " << (A.empty() ? std::string("-inf") : llToStr(degOf(A))) << "]\n";
        std::cout << "B(x) = " << polyToStr(B)
                  << "      [deg B = " << (B.empty() ? std::string("-inf") : llToStr(degOf(B))) << "]\n\n";

        // ---- 带余除法 ----
        std::cout << "【带余除法】\n";
        if (B.empty()) {
            std::cout << "  因为 B(x) 是零多项式，带余除法无意义（除数不能为 0）。\n\n";
        } else {
            Poly Q, R;
            polyDivMod(A, B, Q, R);
            std::cout << "  商     Q(x) = " << polyToStr(Q) << "\n";
            std::cout << "  余式   R(x) = " << polyToStr(R) << "\n";
            std::cout << "  即  A(x) = B(x) * Q(x) + R(x)\n";
            std::cout << "          " << polyToStr(A) << " = (" << polyToStr(B)
                      << ") * (" << polyToStr(Q) << ") + (" << polyToStr(R) << ")\n";
            // 验证
            Poly chk = addPoly(mulPoly(B, Q), R);
            std::cout << "  校验：B*Q + R = " << polyToStr(chk)
                      << (eqPoly(chk, A) ? "   [与 A 相同, 正确]\n" : "   [与 A 不同, 错误!]\n");
            std::cout << "  余式次数 deg R = "
                      << (R.empty() ? std::string("-inf") : llToStr(degOf(R)))
                      << "  <  deg B = " << llToStr(degOf(B)) << "\n\n";
        }

        // ---- 最大公因数 ----
        std::cout << "【最大公因数 gcd(A, B)】\n";
        std::vector<std::string> steps;
        Poly G = polyGcd(A, B, &steps);
        std::cout << "  首一形式 (monic)          gcd = " << polyToStr(G) << "\n";
        std::cout << "  整数本原形式(消分母/提内容) gcd = " << polyToStr(primitiveInteger(G)) << "\n";
        std::cout << "  (域上有理系数下 gcd 相差一个非零常数倍，这里以首一形式为代表)\n";

        if (A.empty() && B.empty())
            std::cout << "  注：A、B 都是零多项式，gcd 规定为 0。\n";

        // ---- 裴蜀等式 u*A + v*B = gcd ----
        std::cout << "\n【裴蜀等式】求 u(x), v(x) 使 u(x)·A(x) + v(x)·B(x) = gcd(A, B)\n";
        if (A.empty() && B.empty()) {
            std::cout << "  因为 A = B = 0，等式成为 0 = 0，u(x)、v(x) 可任取。\n";
        } else {
            Poly U, V, D;
            polyExtGcd(A, B, U, V, D);
            std::cout << "  u(x) = " << polyToStr(U) << "\n";
            std::cout << "  v(x) = " << polyToStr(V) << "\n";
            Poly bez = addPoly(mulPoly(U, A), mulPoly(V, B));
            std::cout << "  校验 u·A + v·B = " << polyToStr(bez)
                      << (eqPoly(bez, G) ? "   [= gcd, 正确]\n" : "   [错误!]\n");
            std::cout << "  注：u(x)、v(x) 不唯一；u + k·(B/d)、v - k·(A/d) 仍是解\n";
            std::cout << "      （k(x) 为任意多项式，d = gcd(A, B)）。\n";
        }

        std::cout << "\n";

        // ---- 是否显示中间过程 ----
        std::cout << "是否显示辗转相除的中间过程? (y/N) > ";
        std::string show;
        if (!std::getline(std::cin, show)) break;
        if (!show.empty() && (show[0] == 'y' || show[0] == 'Y')) {
            std::cout << "\n---------- 辗转相除过程 ----------\n";
            if (steps.empty())
                std::cout << "  (B 为零多项式，无需辗转相除)\n";
            for (size_t i = 0; i < steps.size(); ++i) std::cout << steps[i] << "\n";
            std::cout << "-----------------------------------\n";
        }

        if (g_overflow)
            std::cout << "\n[警告] 中间结果超出 64 位整数范围，数值可能不精确。请注意输入规模。\n";

        std::cout << "\n继续计算下一组? (Y/n) > ";
        std::string again;
        if (!std::getline(std::cin, again)) break;
        if (!again.empty() && (again[0] == 'n' || again[0] == 'N')) break;
        std::cout << "\n";
    }

    std::cout << "\n程序结束。\n";
    return 0;
}