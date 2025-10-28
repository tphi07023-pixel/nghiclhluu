#include <math.h>

#define CTRL_SIGN      (+1.0)
#define VDC_EST        (400.0)
#define M_MAX          (0.95)

__declspec(dllexport)
void simuser(t, delta, in, out)
double* in, * out;
double t, delta;
{
    // States
    static double Iint = 0.0, Ramp = 0.0, Ifil = 0.0;
    static double x1 = 0.0, x2 = 0.0, zpll = 0.0, th = 0.0;   // SOGI & góc PLL
    static double Ipk = 0.0, Ipk_prev = 0.0;

    // --- Auto-cal cho c?m bi?n dòng ---
    static double Gcal = 1.0;                 // h? s? scale ?o (i_meas = Gcal*i_raw)
    static double idr_lp = 0.0, iqr_lp = 0.0; // LPF tách sóng ??ng b?

    // Tunings (vòng dòng)
    const double Ipk_set = 10.0;
    const double f_carrier = 10000.0;
    const double Kp = 0.03;
    const double Ki = 12.0;
    const double Kaw = 120.0;
    const double alpha_I = 0.20;
    const double Ipk_slew = 800.0;
    const double TWO_PI = 6.283185307179586;
    const double Ts = delta;

    // Plant & FF
    const double L_EQ = 20e-3;   // L2
    const double KffL = 0.80;
    const double KffV = 0.90;

    // Inputs
    const double i_raw = in[0];   // dòng th?t t? mô ph?ng (ch?a scale)
    const double v_pcc = in[1];

    // --- SOGI-PLL ---
    const double w0 = TWO_PI * 50.0, k = 1.0;
    x1 += Ts * (-k * w0 * x1 - w0 * w0 * x2 + k * w0 * v_pcc);
    x2 += Ts * (x1);
    const double valpha = x1;
    const double vbeta = w0 * x2;

    const double s = sin(th), c = cos(th);
    const double vd = valpha * s - vbeta * c;     // ~ DC d??ng khi khóa
    const double vq = valpha * c + vbeta * s;     // ~ 0 khi khóa

    // PLL chu?n hóa + anti-windup
    double Vpk = sqrt(valpha * valpha + vbeta * vbeta);
    if (Vpk < 10.0) Vpk = 10.0;
    const double vq_n = vq / Vpk;

    const double zeta = 0.707, w_pll = TWO_PI * 30.0;
    const double Kp_pll = 2.0 * zeta * w_pll;
    const double Ki_pll = w_pll * w_pll;

    double w_cmd = w0 + (Kp_pll * vq_n + zpll);
    double w_hat = w_cmd;
    const double wmin = TWO_PI * 49.0, wmax = TWO_PI * 51.0;
    if (w_hat < wmin) w_hat = wmin;
    if (w_hat > wmax) w_hat = wmax;
    zpll += Ki_pll * Ts * (vq_n + (w_hat - w_cmd) / Kp_pll);

    th += w_hat * Ts;
    if (th >= TWO_PI) th -= TWO_PI;
    if (th < 0.0)     th += TWO_PI;

    // --- Headroom & ramp Ipk ---
    const double Vg_pk_hat = Vpk;
    double Ipk_max = (M_MAX * VDC_EST - Vg_pk_hat) / (w0 * L_EQ);
    if (Ipk_max < 0.0) Ipk_max = 0.0;
    double Ipk_target = Ipk_set;
    if (Ipk_target > 0.9 * Ipk_max) Ipk_target = 0.9 * Ipk_max;

    double dI = Ipk_target - Ipk;
    double lim = Ipk_slew * Ts;
    if (dI > lim) dI = lim;
    if (dI < -lim) dI = -lim;
    Ipk += dI;

    // --- Auto-cal: ??c l??ng biên ?? i_raw & c?p nh?t Gcal ---
    const double w_agc = TWO_PI * 10.0;         // BW ~10 Hz
    idr_lp += Ts * (w_agc * (i_raw * s - idr_lp));
    iqr_lp += Ts * (w_agc * (i_raw * c - iqr_lp));
    double Ipk_raw = 2.0 * sqrt(idr_lp * idr_lp + iqr_lp * iqr_lp) + 1e-6;

    double Gcal_tgt = Ipk / Ipk_raw;            // m?c tiêu ?? i_meas_pk = Ipk
    double mu = 0.05;                           // t?c ?? h?c
    Gcal += mu * (Gcal_tgt - Gcal);
    if (Gcal < 0.3) Gcal = 0.3;
    if (Gcal > 3.0) Gcal = 3.0;

    // Dòng ??a vào vòng dòng
    const double i_meas = Gcal * i_raw;
    Ifil = (1.0 - alpha_I) * Ifil + alpha_I * i_meas;

    // --- Iref & ??o hàm s?ch ---
    const double Iref = Ipk * s;
    const double dIpk = (Ipk - Ipk_prev) / Ts;  Ipk_prev = Ipk;
    const double dIref = Ipk * w_hat * c + dIpk * s;

    // --- PI dòng + feed-forward ---
    const double e = Iref - Ifil;
    const double vff_grid = KffV * (vd / VDC_EST);
    const double vff_Ldi = (KffL * L_EQ / VDC_EST) * dIref;

    const double u_pi = CTRL_SIGN * (Kp * e + Iint);
    const double u = u_pi + vff_grid + vff_Ldi;

    // Saturation & anti-windup
    double m = u;
    if (m > M_MAX) m = M_MAX;
    if (m < -M_MAX) m = -M_MAX;
    int sat_hi = (u > M_MAX), sat_lo = (u < -M_MAX);
    int pushing = ((sat_hi && (CTRL_SIGN * e > 0)) || (sat_lo && (CTRL_SIGN * e < 0)));
    if (pushing) Iint += Kaw * (m - u) * Ts;
    else         Iint += (Ki * e + Kaw * (m - u)) * Ts;

    // PWM
    double Duty = 0.5 * (m + 1.0);
    Ramp += Ts * f_carrier; if (Ramp >= 1.0) Ramp -= 1.0;
    double Sa = (Ramp < Duty) ? 1.0 : 0.0;

    // Outputs
    out[0] = Sa;
    out[1] = Duty;
    out[2] = Iref;
    out[3] = vq;
    out[4] = w_hat / TWO_PI;   // Hz
    out[5] = u;                // tr??c k?p
    out[6] = m;                // sau k?p
    out[7] = e;                // sai s?
    out[8] = Ipk;
    out[9] = Ipk_max;
    out[10] = Ipk_raw;          // biên ?? dòng RAW
    out[11] = Gcal;             // h? s? scale ?ang h?c
    out[12] = Gcal_tgt;         // m?c tiêu Gcal
}
