// ============================================================================
//  poly_ui.cpp  --  一元多项式 带余除法 + 最大公因数  (EasyX 图形界面版)
//
//  依赖：EasyX for MinGW（graphics.h + libeasyxw.a）
//  编译：
//      g++ -std=c++11 -O2 -DUNICODE -D_UNICODE -mwindows poly_ui.cpp -o poly_ui.exe ^
//          -leasyxw -lgdi32 -luser32 -lole32 -luuid -loleaut32 -limm32 -lwinmm
//
//  关于清晰度：
//    * 程序启动即声明「DPI 感知」(SetProcessDPIAware)，避免 Windows 对窗口做
//      2 倍位图拉伸（那是高分屏发虚的根因）。
//    * 所有界面坐标按 96 DPI 的“设计像素”书写，绘制前用 setaspectratio(scale,scale)
//      换算到物理像素，文字由 GDI 在物理分辨率上重新绘制 —— 真·高清，不是拉大。
//    * 字体统一用自建 LOGFONT + CLEARTYPE_QUALITY 渲染，边缘平滑。
// ============================================================================

// 需要 Win7+ 的 API 声明（SetProcessDPIAware 等）
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif
#ifndef WINVER
#define WINVER 0x0601
#endif

#include <graphics.h>
#include <windows.h>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <climits>

// ===========================================================================
//  一、有理数（分数）与多项式运算内核
// ===========================================================================
#if defined(__SIZEOF_INT128__)
typedef __int128 wide_t;
#define HAS_WIDE 1
#else
typedef long long wide_t;
#define HAS_WIDE 0
#endif

static bool g_overflow = false;

static wide_t wabs(wide_t x) { return x < 0 ? -x : x; }
static wide_t wgcd(wide_t a, wide_t b) {
    a = wabs(a); b = wabs(b);
    while (b != 0) { wide_t t = a % b; a = b; b = t; }
    return a;
}
static long long toLL(wide_t x) {
#if HAS_WIDE
    if (x > (wide_t)LLONG_MAX || x < (wide_t)LLONG_MIN) { g_overflow = true; return (long long)(x < 0 ? LLONG_MIN : LLONG_MAX); }
#endif
    return (long long)x;
}

struct Frac {
    long long p;   // 分子
    long long q;   // 分母 (>0)
    Frac() : p(0), q(1) {}
    explicit Frac(long long v) : p(0), q(1) { set((wide_t)v, (wide_t)1); }
    Frac(long long pp, long long qq) : p(0), q(1) { set((wide_t)pp, (wide_t)qq); }
    void set(wide_t pp, wide_t qq) {
        if (qq == 0) { pp = 0; qq = 1; }
        if (qq < 0) { pp = -pp; qq = -qq; }
        if (pp == 0) { p = 0; q = 1; return; }
        wide_t g = wgcd(pp, qq);
        if (g > 1) { pp /= g; qq /= g; }
        p = toLL(pp); q = toLL(qq);
    }
    bool isZero() const { return p == 0; }
};

static Frac addF(const Frac& a, const Frac& b) { Frac r; r.set((wide_t)a.p * b.q + (wide_t)b.p * a.q, (wide_t)a.q * b.q); return r; }
static Frac subF(const Frac& a, const Frac& b) { Frac r; r.set((wide_t)a.p * b.q - (wide_t)b.p * a.q, (wide_t)a.q * b.q); return r; }
static Frac mulF(const Frac& a, const Frac& b) {
    wide_t n1 = a.p, d1 = a.q, n2 = b.p, d2 = b.q;
    wide_t g1 = wgcd(n1, d2); if (g1 > 1) { n1 /= g1; d2 /= g1; }
    wide_t g2 = wgcd(n2, d1); if (g2 > 1) { n2 /= g2; d1 /= g2; }
    Frac r; r.set(n1 * n2, d1 * d2); return r;
}
static Frac divF(const Frac& a, const Frac& b) {
    wide_t num = a.p, den = a.q, bn = b.p, bd = b.q;
    wide_t g1 = wgcd(num, bn); if (g1 > 1) { num /= g1; bn /= g1; }
    wide_t g2 = wgcd(bd, den); if (g2 > 1) { bd /= g2; den /= g2; }
    Frac r; r.set(num * bd, den * bn); return r;
}

typedef std::vector<Frac> Poly;

static void trim(Poly& a) { while (!a.empty() && a.back().isZero()) a.pop_back(); }
static int degOf(const Poly& a) { return (int)a.size() - 1; }

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

static bool polyDivMod(const Poly& a, const Poly& b, Poly& q, Poly& r) {
    if (b.empty()) return false;
    Poly rem = a;
    int db = (int)b.size() - 1, da = (int)a.size() - 1;
    q.assign((size_t)std::max(0, da - db + 1), Frac(0));
    for (int i = da; i >= db; --i) {
        if (i >= (int)rem.size() || rem[i].isZero()) continue;
        Frac factor = divF(rem[i], b[db]);
        int shift = i - db;
        q[shift] = factor;
        for (int j = 0; j <= db; ++j)
            rem[shift + j] = subF(rem[shift + j], mulF(factor, b[j]));
    }
    trim(q); trim(rem); r = rem;
    return true;
}

static Poly monic(const Poly& a) {
    Poly r = a; trim(r);
    if (r.empty()) return r;
    Frac lead = r.back();
    for (size_t i = 0; i < r.size(); ++i) r[i] = divF(r[i], lead);
    return r;
}

static Poly primitiveInteger(const Poly& a) {
    Poly r = a; trim(r);
    if (r.empty()) return r;
    wide_t l = 1;
    for (size_t i = 0; i < r.size(); ++i) { wide_t g = wgcd(l, r[i].q); l = l / (g ? g : 1) * r[i].q; }
    std::vector<wide_t> num(r.size());
    for (size_t i = 0; i < r.size(); ++i) num[i] = (wide_t)r[i].p * (l / r[i].q);
    wide_t g = 0;
    for (size_t i = 0; i < num.size(); ++i) g = wgcd(g, num[i]);
    if (g == 0) g = 1;
    if (num.back() < 0) g = -g;
    Poly out;
    for (size_t i = 0; i < num.size(); ++i) out.push_back(Frac(toLL(num[i] / g)));
    trim(out); return out;
}

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
        Frac mag = neg ? Frac(-f.p, f.q) : f;
        if (s.empty()) { if (neg) s += "-"; }
        else s += (neg ? " - " : " + ");
        if (i == 0) { s += fracToStr(mag); }
        else {
            bool isOne = (mag.p == 1 && mag.q == 1);
            if (!isOne) { if (mag.q == 1) s += llToStr(mag.p); else s += "(" + fracToStr(mag) + ")"; }
            s += "x";
            if (i > 1) { s += "^"; s += llToStr(i); }
        }
    }
    return s.empty() ? "0" : s;
}

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
            os << "      A = " << polyToStr(a) << "\n";
            os << "      B = " << polyToStr(b) << "\n";
            os << "      Q = " << polyToStr(q) << "\n";
            os << "      R = " << polyToStr(r);
            log->push_back(os.str());
        }
        a = b; b = r;
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

// ---- 系数解析（整数 / 分数 / 小数）----
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
static bool parseFrac(const std::string& tk, Frac& out) {
    if (tk.empty()) return false;
    size_t slash = tk.find('/');
    if (slash != std::string::npos) {
        long long a, b;
        if (!parseIntStr(tk.substr(0, slash), a)) return false;
        if (!parseIntStr(tk.substr(slash + 1), b)) return false;
        if (b == 0) return false;
        out = Frac(a, b); return true;
    }
    size_t dot = tk.find('.');
    if (dot != std::string::npos) {
        std::string digits = tk.substr(0, dot) + tk.substr(dot + 1);
        bool neg = false;
        if (!digits.empty() && (digits[0] == '+' || digits[0] == '-')) { neg = (digits[0] == '-'); digits = digits.substr(1); }
        if (digits.empty()) return false;
        long long num = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            if (!isdigit((unsigned char)digits[i])) return false;
            num = num * 10 + (digits[i] - '0');
        }
        long long den = 1;
        for (size_t i = dot + 1; i < tk.size(); ++i) den *= 10;
        if (neg) num = -num;
        out = Frac(num, den); return true;
    }
    long long a;
    if (!parseIntStr(tk, a)) return false;
    out = Frac(a); return true;
}
static bool parsePolyLine(const std::string& line, Poly& out, std::string& err) {
    std::istringstream iss(line);
    std::string tk;
    std::vector<Frac> cs;
    while (iss >> tk) {
        Frac f;
        if (!parseFrac(tk, f)) { err = "无法识别的系数: \"" + tk + "\""; return false; }
        cs.push_back(f);
    }
    std::reverse(cs.begin(), cs.end());
    out = cs; trim(out);
    return true;
}

// ===========================================================================
//  二、字符串 & 高清文字工具
// ===========================================================================
static std::wstring widen(const std::string& s) {
    std::wstring w; w.resize(s.size());
    for (size_t i = 0; i < s.size(); ++i) w[i] = (wchar_t)(unsigned char)s[i];
    return w;
}
static std::string narrow(const std::wstring& s) {
    std::string r;
    for (size_t i = 0; i < s.size(); ++i) { wchar_t c = s[i]; r += (c < 128) ? (char)c : '?'; }
    return r;
}
static std::vector<std::wstring> toLines(const std::wstring& s) {
    std::vector<std::wstring> v; std::wstring cur;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == L'\n') { v.push_back(cur); cur.clear(); }
        else if (s[i] != L'\r') cur += s[i];
    }
    v.push_back(cur);
    return v;
}

// 统一字体入口：CLEARTYPE 抗锯齿，负 lfHeight 表示字符高度（em）
static void setFont(int px, bool bold = false, const wchar_t* face = L"微软雅黑") {
    LOGFONTW lf;
    ZeroMemory(&lf, sizeof(lf));
    lf.lfHeight = -px;
    lf.lfWidth = 0;
    lf.lfWeight = bold ? FW_BOLD : FW_NORMAL;
    lf.lfItalic = FALSE;
    lf.lfUnderline = FALSE;
    lf.lfStrikeOut = FALSE;
    lf.lfCharSet = DEFAULT_CHARSET;
    lf.lfOutPrecision = OUT_TT_PRECIS;
    lf.lfClipPrecision = CLIP_DEFAULT_PRECIS;
    lf.lfQuality = CLEARTYPE_QUALITY;
    lf.lfPitchAndFamily = DEFAULT_PITCH | FF_DONTCARE;
    lstrcpynW(lf.lfFaceName, face, LF_FACESIZE);
    settextstyle(&lf);
}

// ===========================================================================
//  三、界面（设计坐标 = 96 DPI 下的像素）
// ===========================================================================
static const int WIN_W = 1060, WIN_H = 790, MARGIN = 22;

static const COLORREF C_BG      = RGB(240, 243, 248);
static const COLORREF C_PANEL   = RGB(255, 255, 255);
static const COLORREF C_BORDER  = RGB(198, 206, 216);
static const COLORREF C_TEXT    = RGB(30, 38, 50);
static const COLORREF C_GRAY    = RGB(122, 132, 148);
static const COLORREF C_ACCENT  = RGB(0, 112, 206);
static const COLORREF C_ACCENT2 = RGB(12, 140, 100);
static const COLORREF C_HOVER   = RGB(224, 238, 252);

enum { BTN_A, BTN_B, BTN_SWAP, BTN_CLEAR, BTN_DEMO, BTN_CALC, BTN_STEPS, BTN_EXIT };
struct Button { int id, x, y, w, h; std::wstring label; };

static Poly gA, gB;
static std::wstring gRawA, gRawB, gStatus;
static std::vector<std::wstring> gLines, gStepLines;
static bool gShowSteps = false, gHasResult = false, gImageMode = false;
static int  gScroll = 0, gContentH = 0;
static int  gMouseX = -1, gMouseY = -1;
static int  gDemoIdx = 0;

static const int BAR_A_T = 82,  BAR_B_T = 144, BAR_H = 54;
static const int BTN_Y = 212, BTN_H = 46;
static const int PANEL_T = 274, PANEL_B = WIN_H - MARGIN;
static const int CONTENT_T = 316, CONTENT_B = PANEL_B - 42;
static const int LH = 25;

static Button mkBtn(int id, int x, int y, int w, int h, const std::wstring& s) {
    Button b; b.id = id; b.x = x; b.y = y; b.w = w; b.h = h; b.label = s; return b;
}
static void buildButtons(std::vector<Button>& bs) {
    bs.clear();
    bs.push_back(mkBtn(BTN_A,      22, BTN_Y, 132, BTN_H, L"输入 A(x)"));
    bs.push_back(mkBtn(BTN_B,     164, BTN_Y, 132, BTN_H, L"输入 B(x)"));
    bs.push_back(mkBtn(BTN_SWAP,  306, BTN_Y, 112, BTN_H, L"交换 A/B"));
    bs.push_back(mkBtn(BTN_CLEAR, 428, BTN_Y,  92, BTN_H, L"清空"));
    bs.push_back(mkBtn(BTN_DEMO,  530, BTN_Y,  92, BTN_H, L"示例"));
    bs.push_back(mkBtn(BTN_CALC,  632, BTN_Y, 132, BTN_H, L"计算"));
    bs.push_back(mkBtn(BTN_STEPS, 774, BTN_Y, 152, BTN_H, gShowSteps ? L"隐藏过程" : L"显示过程"));
    bs.push_back(mkBtn(BTN_EXIT,  936, BTN_Y, 102, BTN_H, L"退出"));
}
static bool inBtn(const Button& b, int x, int y) {
    return x >= b.x && x < b.x + b.w && y >= b.y && y < b.y + b.h;
}
static void drawPanel(int l, int t, int r, int b, COLORREF fill) {
    setfillcolor(fill); setlinecolor(C_BORDER);
    fillrectangle(l, t, r, b);
}
static void drawButton(const Button& b, bool hover) {
    bool primary = (b.id == BTN_CALC);
    bool active  = (b.id == BTN_STEPS && gShowSteps);
    COLORREF fill = C_PANEL, border = C_BORDER, txt = C_TEXT;
    if (primary) { fill = C_ACCENT; border = C_ACCENT; txt = RGB(255, 255, 255); }
    if (active)  { fill = C_HOVER;  border = C_ACCENT; txt = C_ACCENT; }
    if (hover)   { fill = primary ? RGB(22, 136, 232) : C_HOVER; border = C_ACCENT; if (!primary) txt = C_ACCENT; }

    setfillcolor(fill); setlinecolor(border);
    fillroundrect(b.x, b.y, b.x + b.w, b.y + b.h, 12, 12);

    setFont(18);
    settextcolor(txt);
    int tw = textwidth(b.label.c_str());
    int th = textheight(b.label.c_str());
    outtextxy(b.x + (b.w - tw) / 2, b.y + (b.h - th) / 2, b.label.c_str());
}
static std::wstring fitText(const std::wstring& s, int maxW) {
    if (textwidth(s.c_str()) <= maxW) return s;
    std::wstring t = s;
    while (!t.empty() && textwidth((t + L"...").c_str()) > maxW) t.erase(t.size() - 1);
    return t + L"...";
}
static void drawPolyBar(int t, const std::wstring& name, const Poly& p, bool entered, COLORREF accent) {
    const int l = MARGIN, r = WIN_W - MARGIN, b = t + BAR_H;
    drawPanel(l, t, r, b, C_PANEL);
    setfillcolor(accent);
    solidrectangle(l, t, l + 6, b);

    setFont(22, true);
    settextcolor(accent);
    outtextxy(l + 20, t + 15, name.c_str());
    int nameW = textwidth(name.c_str());

    std::wstring shown = entered ? widen(polyToStr(p)) : std::wstring(L"(未输入，点击「输入」按钮)");
    setFont(24);
    settextcolor(entered ? C_TEXT : C_GRAY);
    int x = l + 20 + nameW + 16;
    outtextxy(x, t + 13, fitText(shown, r - x - 130).c_str());

    if (entered) {
        std::wstring d = p.empty() ? std::wstring(L"deg = -inf") : (L"deg = " + widen(llToStr(degOf(p))));
        setFont(15);
        settextcolor(C_GRAY);
        int dw = textwidth(d.c_str());
        outtextxy(r - 18 - dw, t + 20, d.c_str());
    }
}

static std::vector<std::wstring> allLines() {
    std::vector<std::wstring> v = gLines;
    if (gShowSteps && !gStepLines.empty()) {
        v.push_back(L"");
        v.push_back(L"-------- 辗转相除过程 --------");
        for (size_t i = 0; i < gStepLines.size(); ++i) v.push_back(gStepLines[i]);
    }
    return v;
}

static void drawAll() {
    std::vector<Button> buttons;
    buildButtons(buttons);

    setbkcolor(C_BG);
    cleardevice();
    setbkmode(TRANSPARENT);

    // 标题
    setFont(30, true);
    settextcolor(RGB(24, 34, 50));
    outtextxy(MARGIN, 14, L"一元多项式计算器");
    setFont(15);
    settextcolor(C_GRAY);
    outtextxy(MARGIN + 2, 54, L"带余除法 A = B·Q + R     最大公因数 gcd(A, B)     裴蜀等式 u·f + v·g = gcd     系数支持分数/小数");

    drawPolyBar(BAR_A_T, L"A(x) =", gA, !gRawA.empty(), C_ACCENT);
    drawPolyBar(BAR_B_T, L"B(x) =", gB, !gRawB.empty(), C_ACCENT2);

    for (size_t i = 0; i < buttons.size(); ++i)
        drawButton(buttons[i], inBtn(buttons[i], gMouseX, gMouseY));

    // 结果面板
    drawPanel(MARGIN, PANEL_T, WIN_W - MARGIN, PANEL_B, C_PANEL);
    setFont(18, true);
    settextcolor(RGB(24, 34, 50));
    outtextxy(MARGIN + 18, PANEL_T + 10, L"计算结果");
    setlinecolor(RGB(228, 233, 240));
    line(MARGIN + 18, PANEL_T + 38, WIN_W - MARGIN - 18, PANEL_T + 38);

    std::vector<std::wstring> lines = allLines();
    int viewH = CONTENT_B - CONTENT_T;
    gContentH = (int)lines.size() * LH;
    int maxScroll = std::max(0, gContentH - viewH);
    if (gScroll > maxScroll) gScroll = maxScroll;
    if (gScroll < 0) gScroll = 0;

    setFont(19);
    if (lines.empty()) {
        settextcolor(C_GRAY);
        outtextxy(MARGIN + 22, CONTENT_T + 4, L"点击「输入 A(x)」「输入 B(x)」录入多项式，再点「计算」。");
        outtextxy(MARGIN + 22, CONTENT_T + 4 + LH, L"也可以直接点「示例」看一个演示。");
    } else {
        for (size_t i = 0; i < lines.size(); ++i) {
            int y = CONTENT_T + (int)i * LH - gScroll;
            if (y < CONTENT_T - LH / 2 || y + LH > CONTENT_B) continue;
            const std::wstring& s = lines[i];
            settextcolor(C_TEXT);
            if (!s.empty() && s[0] == L'[') settextcolor(C_ACCENT);
            if (s.size() > 2 && s[0] == L'-' && s[1] == L'-') settextcolor(C_GRAY);
            if (s.size() > 4 && s.compare(0, 4, L"  商 ") == 0) settextcolor(C_ACCENT);
            if (s.size() > 5 && s.compare(0, 5, L"  余式") == 0) settextcolor(C_ACCENT2);
            outtextxy(MARGIN + 22, y, fitText(s, WIN_W - 2 * MARGIN - 66).c_str());
        }
        if (maxScroll > 0) {
            int barH = std::max(34, viewH * viewH / gContentH);
            int barT = CONTENT_T + (viewH - barH) * gScroll / maxScroll;
            setfillcolor(RGB(222, 228, 238));
            solidrectangle(WIN_W - MARGIN - 13, CONTENT_T, WIN_W - MARGIN - 7, CONTENT_B);
            setfillcolor(C_ACCENT);
            solidrectangle(WIN_W - MARGIN - 13, barT, WIN_W - MARGIN - 7, barT + barH);
        }
    }

    // 状态栏
    setlinecolor(RGB(228, 233, 240));
    line(MARGIN + 18, PANEL_B - 30, WIN_W - MARGIN - 18, PANEL_B - 30);
    setFont(15);
    settextcolor(C_GRAY);
    outtextxy(MARGIN + 18, PANEL_B - 25, fitText(gStatus, WIN_W - 2 * MARGIN - 44).c_str());
}

static void render() {
    if (gImageMode) { drawAll(); return; }
    BeginBatchDraw();
    drawAll();
    EndBatchDraw();
}

// ===========================================================================
//  四、动作
// ===========================================================================
static void compute() {
    std::vector<std::wstring> out;
    gStepLines.clear();


    if (gB.empty()) {
        out.push_back(L"[带余除法]");
        out.push_back(L"  B(x) 是零多项式，不能作除数，带余除法无意义。");
    } else {
        Poly Q, R;
        polyDivMod(gA, gB, Q, R);
        out.push_back(L"[带余除法]  A = B·Q + R ， deg R < deg B");
        out.push_back(L"  商   Q(x) = " + widen(polyToStr(Q)));
        out.push_back(L"  余式 R(x) = " + widen(polyToStr(R)));
        out.push_back(L"  即： " + widen(polyToStr(gA)) + L" = (" + widen(polyToStr(gB)) + L") · ("
                      + widen(polyToStr(Q)) + L") + (" + widen(polyToStr(R)) + L")");
        Poly chk = addPoly(mulPoly(gB, Q), R);
        out.push_back(std::wstring(L"  校验 B·Q + R = ") + widen(polyToStr(chk))
                      + (eqPoly(chk, gA) ? L"     [正确]" : L"     [错误!]"));
    }
    out.push_back(L"");

    std::vector<std::string> log;
    Poly G = polyGcd(gA, gB, &log);
    out.push_back(L"[最大公因数]");
    out.push_back(L"  首一形式 (monic)               gcd = " + widen(polyToStr(G)));
    out.push_back(L"  整数本原形式 (消分母、提内容)    gcd = " + widen(polyToStr(primitiveInteger(G))));
    out.push_back(L"  注：有理系数下 gcd 相差一个非零常数倍，此处以首一形式为代表。");
    if (gA.empty() && gB.empty()) out.push_back(L"  注：A、B 都是零多项式，规定 gcd = 0。");
    out.push_back(L"");

    out.push_back(L"[裴蜀等式]  u(x)·A(x) + v(x)·B(x) = gcd(A, B)");
    if (gA.empty() && gB.empty()) {
        out.push_back(L"  A = B = 0，等式成为 0 = 0，u(x)、v(x) 可任取。");
    } else {
        Poly U, V, D;
        polyExtGcd(gA, gB, U, V, D);
        out.push_back(L"  u(x) = " + widen(polyToStr(U)));
        out.push_back(L"  v(x) = " + widen(polyToStr(V)));
        Poly bez = addPoly(mulPoly(U, gA), mulPoly(V, gB));
        out.push_back(std::wstring(L"  校验 u·A + v·B = ") + widen(polyToStr(bez))
                      + (eqPoly(bez, G) ? L"     [= gcd, 正确]" : L"     [错误!]"));
        out.push_back(L"  注：u、v 不唯一，u + k·(B/d)、v - k·(A/d) 也是解（k 为任意多项式）。");
    }
    out.push_back(L"");

    if (g_overflow) out.push_back(L"[警告] 中间结果超出 64 位整数范围，数值可能不精确。");

    gLines = out;
    for (size_t i = 0; i < log.size(); ++i) {
        std::vector<std::wstring> ls = toLines(widen(log[i]));
        for (size_t j = 0; j < ls.size(); ++j) gStepLines.push_back(ls[j]);
    }
    gHasResult = true;
    gScroll = 0;
    gStatus = L"计算完成。可点「显示过程」查看辗转相除的每一步。";
}

static void inputPoly(bool isA) {
    wchar_t buf[4096];
    const wchar_t* cur = isA ? gRawA.c_str() : gRawB.c_str();
    lstrcpynW(buf, cur, 4096);
    if (buf[0] == 0) lstrcpynW(buf, isA ? L"1 0 0 -1" : L"1 -1", 4096);

    const wchar_t* title = isA ? L"输入 A(x) 的系数" : L"输入 B(x) 的系数";
    const wchar_t* prompt = L"按【从高次到低次】输入全部系数，空格分隔。\n例：1 -2 3/4 5  表示  x^3 - 2x^2 + (3/4)x + 5\n支持整数/分数(3/4)/小数(0.75)；输入 0 表示零多项式。";

    if (!InputBox(buf, 4096, prompt, title, buf, 640, 230, true)) return;

    std::string err;
    Poly P;
    if (!parsePolyLine(narrow(buf), P, err)) {
        std::wstring m = L"输入有误：" + widen(err) + L"\n\n请检查系数格式（如：1 -2 3/4 5）。";
        MessageBoxW(GetHWnd(), m.c_str(), L"格式错误", MB_OK | MB_ICONWARNING);
        return;
    }
    if (isA) { gA = P; gRawA = buf; gStatus = L"A(x) 已更新。"; }
    else     { gB = P; gRawB = buf; gStatus = L"B(x) 已更新。"; }
    gHasResult = false; gLines.clear(); gStepLines.clear(); gScroll = 0;
}

static void loadDemo() {
    struct Demo { const wchar_t* a; const wchar_t* b; const wchar_t* tip; };
    static const Demo demos[] = {
        { L"1 0 0 0 -1", L"1 -2 -1 2", L"示例 1：x^4-1 与 x^3-2x^2-x+2，gcd = x^2-1" },
        { L"2 3 1",      L"1 1/2",     L"示例 2：含分数系数的整除" },
        { L"6 -5 -2 1",  L"3 1",       L"示例 3：带余除法出现分数商" }
    };
    const int n = (int)(sizeof(demos) / sizeof(demos[0]));
    const Demo& d = demos[gDemoIdx % n];
    gDemoIdx = (gDemoIdx + 1) % n;

    std::string err;
    Poly Pa, Pb;
    parsePolyLine(narrow(d.a), Pa, err);
    parsePolyLine(narrow(d.b), Pb, err);
    gA = Pa; gB = Pb; gRawA = d.a; gRawB = d.b;
    gHasResult = false; gLines.clear(); gStepLines.clear(); gScroll = 0;
    gStatus = d.tip;
}

static void handleClick(int x, int y) {
    std::vector<Button> bs;
    buildButtons(bs);
    for (size_t i = 0; i < bs.size(); ++i) {
        if (!inBtn(bs[i], x, y)) continue;
        switch (bs[i].id) {
        case BTN_A:     inputPoly(true);  break;
        case BTN_B:     inputPoly(false); break;
        case BTN_SWAP:
            std::swap(gA, gB); std::swap(gRawA, gRawB);
            gStatus = L"已交换 A(x) 与 B(x)。";
            gHasResult = false; gLines.clear(); gStepLines.clear(); gScroll = 0;
            break;
        case BTN_CLEAR:
            gA.clear(); gB.clear(); gRawA.clear(); gRawB.clear();
            gLines.clear(); gStepLines.clear(); gHasResult = false; gShowSteps = false; gScroll = 0;
            gStatus = L"已清空。"; break;
        case BTN_DEMO:  loadDemo(); break;
        case BTN_CALC:
            if (gRawA.empty() && gRawB.empty()) { gStatus = L"请先输入 A(x) 和 B(x)。"; break; }
            compute(); break;
        case BTN_STEPS:
            if (!gHasResult) { gStatus = L"请先点「计算」。"; break; }
            gShowSteps = !gShowSteps; gScroll = 0;
            gStatus = gShowSteps ? L"已展开辗转相除过程。" : L"已收起辗转相除过程。";
            break;
        case BTN_EXIT:  closegraph(); exit(0); break;
        }
        return;
    }
}

// ===========================================================================
//  五、DPI 缩放 & 主循环
// ===========================================================================
static float gScale = 1.0f;            // 物理像素 / 设计像素
static int   SC(int v) { return (int)(v * gScale + 0.5f); }

static void enableDpiAwareness() {
    // 这些 API 在本机旧版 MinGW 头文件里没有声明，一律用 GetProcAddress 动态获取。
    HMODULE user32 = LoadLibraryW(L"user32.dll");
    if (!user32) return;

    // 优先 Per-Monitor V2（Win10 1703+），DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2 == (void*)-4
    typedef BOOL (WINAPI *SetCtxFn)(void*);
    SetCtxFn setCtx = (SetCtxFn)(void*)GetProcAddress(user32, "SetProcessDpiAwarenessContext");
    if (setCtx && setCtx((void*)-4)) { FreeLibrary(user32); return; }

    // 其次 System DPI Aware（Vista+）
    typedef BOOL (WINAPI *SetAwareFn)(void);
    SetAwareFn setAware = (SetAwareFn)(void*)GetProcAddress(user32, "SetProcessDPIAware");
    if (setAware) setAware();

    FreeLibrary(user32);
}
static void computeScale() {
    HDC dc = GetDC(NULL);
    int dpi = dc ? GetDeviceCaps(dc, LOGPIXELSX) : 96;
    if (dc) ReleaseDC(NULL, dc);
    if (dpi < 96) dpi = 96;
    gScale = dpi / 96.0f;
}

static int runApp() {
    initgraph(SC(WIN_W), SC(WIN_H));
    setaspectratio(gScale, gScale);      // 设计坐标 -> 物理像素（文字按真实分辨率重绘）
    HWND hwnd = GetHWnd();
    SetWindowTextW(hwnd, L"一元多项式计算器  -  带余除法 & 最大公因数 & 裴蜀等式");

    // 设置窗口 / 任务栏 / Alt-Tab 图标（图标资源 id = 1，见 app.rc）
    HINSTANCE hInst = GetModuleHandleW(NULL);
    HICON hBig = (HICON)LoadImageW(hInst, MAKEINTRESOURCEW(1), IMAGE_ICON,
                                   GetSystemMetrics(SM_CXICON), GetSystemMetrics(SM_CYICON), 0);
    HICON hSmall = (HICON)LoadImageW(hInst, MAKEINTRESOURCEW(1), IMAGE_ICON,
                                     GetSystemMetrics(SM_CXSMICON), GetSystemMetrics(SM_CYSMICON), 0);
    if (hBig)   SendMessageW(hwnd, WM_SETICON, ICON_BIG,   (LPARAM)hBig);
    if (hSmall) SendMessageW(hwnd, WM_SETICON, ICON_SMALL, (LPARAM)hSmall);
    if (hBig)   SetClassLongPtrW(hwnd, GCLP_HICON,   (LONG_PTR)hBig);
    if (hSmall) SetClassLongPtrW(hwnd, GCLP_HICONSM, (LONG_PTR)hSmall);

    LONG st = GetWindowLong(hwnd, GWL_STYLE);
    st &= ~(WS_MAXIMIZEBOX | WS_THICKFRAME);
    SetWindowLong(hwnd, GWL_STYLE, st);
    SetWindowPos(hwnd, NULL, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOZORDER | SWP_FRAMECHANGED);
    SetWindowLong(hwnd, GWL_EXSTYLE, GetWindowLong(hwnd, GWL_EXSTYLE) & ~WS_EX_CLIENTEDGE);

    gStatus = L"就绪：点击「输入 A(x)」「输入 B(x)」录入多项式，或点「示例」。";
    gShowSteps = false;

    bool running = true;
    render();
    while (running) {
        ExMessage msg;
        bool dirty = false;
        // 鼠标坐标是物理像素，换算回设计坐标再命中测试
        while (peekmessage(&msg, EX_MOUSE | EX_KEY | EX_WINDOW)) {
            if (msg.message == WM_LBUTTONDOWN) {
                handleClick((int)(msg.x / gScale), (int)(msg.y / gScale));
                dirty = true;
            } else if (msg.message == WM_MOUSEMOVE) {
                if (msg.x != gMouseX || msg.y != gMouseY) { gMouseX = msg.x; gMouseY = msg.y; dirty = true; }
            } else if (msg.message == WM_MOUSEWHEEL) {
                gScroll -= (int)msg.wheel / 120 * 44;
                dirty = true;
            } else if (msg.message == WM_KEYDOWN) {
                if (msg.vkcode == VK_ESCAPE) running = false;
                else if (msg.vkcode == VK_UP)   { gScroll -= 44; dirty = true; }
                else if (msg.vkcode == VK_DOWN) { gScroll += 44; dirty = true; }
            } else if (msg.message == WM_CLOSE) {
                running = false;
            }
        }
        if (!IsWindow(hwnd)) break;
        if (dirty) render();
        Sleep(10);
    }
    closegraph();
    return 0;
}

static int shotMode(const wchar_t* path, float s) {
    loadDemo(); compute(); gShowSteps = false;
    initgraph(SC(WIN_W), SC(WIN_H));
    IMAGE img((int)(WIN_W * s + 0.5f), (int)(WIN_H * s + 0.5f));
    SetWorkingImage(&img);
    setaspectratio(s, s);                // s=1 出设计稿，s=gScale 出“真实窗口”像素
    gImageMode = true;
    gScale = s;
    drawAll();
    gImageMode = false;
    SetWorkingImage(NULL);
    saveimage(path, &img);
    closegraph();
    return 0;
}

int main(int argc, char** argv) {
    enableDpiAwareness();
    computeScale();
    if (argc >= 3 && std::string(argv[1]) == "--shot") {
        std::string p = argv[2];
        std::wstring wp(p.begin(), p.end());
        float s = (argc >= 4) ? (float)atof(argv[3]) : 1.0f;
        if (s <= 0.1f) s = 1.0f;
        return shotMode(wp.c_str(), s);
    }
    return runApp();
}