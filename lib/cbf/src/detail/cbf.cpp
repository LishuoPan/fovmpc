#include <vector>
#include <cbf/detail/cbf.h>


// Default: alpha = gamma * b
GiNaC::ex defaultAlpha(GiNaC::ex b, double gamma = 1.0)
{
    return gamma * b;
}

GiNaC::ex myAlpha(GiNaC::ex myh, double mygamma)
{
    return mygamma * GiNaC::pow(myh, 5);
}

GiNaC::ex fifthAlpha(GiNaC::ex myh, double mygamma)
{
    return mygamma * GiNaC::pow(myh, 5);
}

namespace cbf
{   
    // template <typename T>
    FovCBF::FovCBF(double fov, double safety_dist, double max_dist, double vmax) : 
    fov(fov), Ds(safety_dist), Rs(max_dist), vmax(vmax), vmin(-vmax), px("px"), py("py"), th("th"), vx("vx"), vy("vy"), w("w"), xt("xt"), yt("yt")
    {
        STATE_VARS = 6;
        CONTROL_VARS = 3;
        gamma = 0.1;

        // GiNaC::symbol px("px"), py("py"), th("th"), vx("vx"), vy("vy"), w("w"), xt("xt"), yt("yt"); //, Ds("Ds"), fov("fov"), Rs("Rs");


        A = {{0, 0, 0, 1, 0, 0},
            {0, 0, 0, 0, 1, 0},
            {0, 0, 0, 0, 0, 1},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0}};

        B = {{0, 0, 0},
            {0, 0, 0},
            {0, 0, 0},
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}};

        state = {{px, py, th, vx, vy, w}};
        state = state.transpose();
        x_target = {{xt, yt}};
        x_target = x_target.transpose();

        f = A.mul(state);
        g = B;

        // Set default alpha
        alpha = fifthAlpha;

        // ------- Safety constraint ------- //
        auto res = initSafetyCBF();
        Ac_safe = res.first;
        Bc_safe = res.second;

        auto res2 = initBorder1CBF();
        Ac_lb = res2.first;
        Bc_lb = res2.second;

        auto res3 = initBorder2CBF();
        Ac_rb = res3.first;
        Bc_rb = res3.second;

        auto res4 = initRangeCBF();
        Ac_range = res4.first;
        Bc_range = res4.second;


        GiNaC::ex bv1 = -state[3] + vmax;
        res = initVelCBF(bv1);
        Ac_v1_max = res.first;
        Bc_v1_max = res.second;

        GiNaC::ex bv2 = -state[4] + vmax;
        res = initVelCBF(bv2);
        Ac_v2_max = res.first;
        Bc_v2_max = res.second;


        GiNaC::ex bv3 = -state[5] + vmax;
        res = initVelCBF(bv3);
        Ac_v3_max = res.first;
        Bc_v3_max = res.second;

        GiNaC::ex bv4 = state[3] - vmin;
        res = initVelCBF(bv4);
        Ac_v1_min = res.first;
        Bc_v1_min = res.second;

        GiNaC::ex bv5 = state[4] - vmin;
        res = initVelCBF(bv5);
        Ac_v2_min = res.first;
        Bc_v2_min = res.second;

        GiNaC::ex bv6 = state[5] - vmin;
        res = initVelCBF(bv6);
        Ac_v3_min = res.first;
        Bc_v3_min = res.second;

    }

    FovCBF::~FovCBF()
    {
        std::cout << "Closing CBF ..." << std::endl;
    }

    
    std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initSafetyCBF()
    {
        GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});
        GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)},
                           {-GiNaC::sin(th), GiNaC::cos(th)}};
        GiNaC::matrix xt_rel = R.mul(d2);
        GiNaC::ex norm2 = GiNaC::pow(xt_rel(0,0), 2) + GiNaC::pow(xt_rel(1,0), 2);

        GiNaC::ex b1 = norm2 - GiNaC::pow(Ds, 2);

        GiNaC::matrix grad_b1 = GiNaC::matrix(STATE_VARS, 1);
        grad_b1(0, 0) = GiNaC::diff(b1, px);
        grad_b1(1, 0) = GiNaC::diff(b1, py);
        grad_b1(2, 0) = GiNaC::diff(b1, th);
        grad_b1(3, 0) = GiNaC::diff(b1, vx);
        grad_b1(4, 0) = GiNaC::diff(b1, vy);
        grad_b1(5, 0) = GiNaC::diff(b1, w);

        GiNaC::ex lfb1 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb1 = lfb1 + grad_b1(i, 0) * f(i, 0);
        }
        
        GiNaC::matrix grad2_b1 = GiNaC::matrix(STATE_VARS, 1);
        grad2_b1(0, 0) = GiNaC::diff(lfb1, px);
        grad2_b1(1, 0) = GiNaC::diff(lfb1, py);
        grad2_b1(2, 0) = GiNaC::diff(lfb1, th);
        grad2_b1(3, 0) = GiNaC::diff(lfb1, vx);
        grad2_b1(4, 0) = GiNaC::diff(lfb1, vy);
        grad2_b1(5, 0) = GiNaC::diff(lfb1, w);
        GiNaC::ex lf2b1 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lf2b1 = lf2b1 + grad2_b1(i, 0) * f(i, 0);
        }

        GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
        GiNaC::ex alpha_b = alpha(b1, gamma);
        grad_bc(0, 0) = GiNaC::diff(alpha_b, px);
        grad_bc(1, 0) = GiNaC::diff(alpha_b, py); 
        grad_bc(2, 0) = GiNaC::diff(alpha_b, th); 
        grad_bc(3, 0) = GiNaC::diff(alpha_b, vx); 
        grad_bc(4, 0) = GiNaC::diff(alpha_b, vy); 
        grad_bc(5, 0) = GiNaC::diff(alpha_b, w);
        GiNaC::ex lfb_c = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
        }

        GiNaC::matrix Ac1 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac1j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac1j += grad2_b1(i, 0) * g(i, j);
            }
            Ac1(0, j) = Ac1j;
        }

        /*
        GiNaC::ex B1 = lf2b1;
        GiNaC::ex B2 = gamma * 3 * GiNaC::pow(b1, 2) * lfb1;
        GiNaC::ex B3 = gamma * GiNaC::pow(lfb1 + gamma*GiNaC::pow(b1, 3), 3);
        GiNaC::ex B4 = 0;
        GiNaC::ex B5 = 0;

        GiNaC::ex Bc1 = B1 + B2 + B3 + B4 + B5;
        */

        GiNaC::ex B1 = lf2b1;
        GiNaC::ex B2 = lfb_c;
        GiNaC::ex B3 = alpha(lfb1 + alpha_b, gamma);
        GiNaC::ex Bc1 = B1 + B2 + B3;
        
        return std::make_pair(Ac1, Bc1);
    }

    std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initBorder1CBF()
    {
        GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});
        GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)},
                            {-GiNaC::sin(th), GiNaC::cos(th)}};
        GiNaC::matrix xt_rel = R.mul(d2);

        GiNaC::ex b2 = GiNaC::tan(fov/2)*xt_rel(0,0) + xt_rel(1,0);
        GiNaC::matrix grad_b2 = GiNaC::matrix(STATE_VARS, 1);
        grad_b2(0, 0) = GiNaC::diff(b2, px);
        grad_b2(1, 0) = GiNaC::diff(b2, py);
        grad_b2(2, 0) = GiNaC::diff(b2, th);
        grad_b2(3, 0) = GiNaC::diff(b2, vx);
        grad_b2(4, 0) = GiNaC::diff(b2, vy);
        grad_b2(5, 0) = GiNaC::diff(b2, w);

        GiNaC::ex lfb2 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb2 = lfb2 + grad_b2(i, 0) * f(i, 0);
        }

        GiNaC::matrix grad2_b2 = GiNaC::matrix(STATE_VARS, 1);
        grad2_b2(0, 0) = GiNaC::diff(lfb2, px);
        grad2_b2(1, 0) = GiNaC::diff(lfb2, py);
        grad2_b2(2, 0) = GiNaC::diff(lfb2, th);
        grad2_b2(3, 0) = GiNaC::diff(lfb2, vx);
        grad2_b2(4, 0) = GiNaC::diff(lfb2, vy);
        grad2_b2(5, 0) = GiNaC::diff(lfb2, w);
        GiNaC::ex lf2b2 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lf2b2 = lf2b2 + grad2_b2(i, 0) * f(i, 0);
        }

        GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
        GiNaC::ex alpha_b2 = alpha(b2, gamma);
        grad_bc(0, 0) = GiNaC::diff(alpha_b2, px);
        grad_bc(1, 0) = GiNaC::diff(alpha_b2, py); 
        grad_bc(2, 0) = GiNaC::diff(alpha_b2, th); 
        grad_bc(3, 0) = GiNaC::diff(alpha_b2, vx); 
        grad_bc(4, 0) = GiNaC::diff(alpha_b2, vy); 
        grad_bc(5, 0) = GiNaC::diff(alpha_b2, w);
        GiNaC::ex lfb_c = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
        }

        GiNaC::matrix Ac2 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac2j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac2j += grad2_b2(i, 0) * g(i, j);
            }
            Ac2(0, j) = Ac2j;
        }


        GiNaC::ex B1 = lf2b2;
        GiNaC::ex B2 = lfb_c;
        GiNaC::ex B3 = alpha(lfb2 + alpha_b2, gamma);
        GiNaC::ex Bc2 = B1 + B2 + B3;
        
        return std::make_pair(Ac2, Bc2);
    }

    std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initBorder2CBF()
    {
        GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});
        GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)},
                            {-GiNaC::sin(th), GiNaC::cos(th)}};
        GiNaC::matrix xt_rel = R.mul(d2);
        GiNaC::ex b3 = GiNaC::tan(fov/2)*xt_rel(0,0) - xt_rel(1,0);
        GiNaC::matrix grad_b3 = GiNaC::matrix(STATE_VARS, 1);
        grad_b3(0, 0) = GiNaC::diff(b3, px);
        grad_b3(1, 0) = GiNaC::diff(b3, py);
        grad_b3(2, 0) = GiNaC::diff(b3, th);
        grad_b3(3, 0) = GiNaC::diff(b3, vx);
        grad_b3(4, 0) = GiNaC::diff(b3, vy);
        grad_b3(5, 0) = GiNaC::diff(b3, w);

        GiNaC::ex lfb3 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb3 = lfb3 + grad_b3(i, 0) * f(i, 0);
        }

        GiNaC::matrix grad2_b3 = GiNaC::matrix(STATE_VARS, 1);
        grad2_b3(0, 0) = GiNaC::diff(lfb3, px);
        grad2_b3(1, 0) = GiNaC::diff(lfb3, py);
        grad2_b3(2, 0) = GiNaC::diff(lfb3, th);
        grad2_b3(3, 0) = GiNaC::diff(lfb3, vx);
        grad2_b3(4, 0) = GiNaC::diff(lfb3, vy);
        grad2_b3(5, 0) = GiNaC::diff(lfb3, w);
        GiNaC::ex lf2b3 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lf2b3 = lf2b3 + grad2_b3(i, 0) * f(i, 0);
        }

        GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
        GiNaC::ex alpha_b3 = alpha(b3, gamma);
        grad_bc(0, 0) = GiNaC::diff(alpha_b3, px);
        grad_bc(1, 0) = GiNaC::diff(alpha_b3, py); 
        grad_bc(2, 0) = GiNaC::diff(alpha_b3, th); 
        grad_bc(3, 0) = GiNaC::diff(alpha_b3, vx); 
        grad_bc(4, 0) = GiNaC::diff(alpha_b3, vy); 
        grad_bc(5, 0) = GiNaC::diff(alpha_b3, w);
        GiNaC::ex lfb_c = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
        }

        GiNaC::matrix Ac3 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac3j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac3j += grad2_b3(i, 0) * g(i, j);
            }
            Ac3(0, j) = Ac3j;
        }


        GiNaC::ex B1 = lf2b3;
        GiNaC::ex B2 = lfb_c;
        GiNaC::ex B3 = alpha(lfb3 + alpha_b3, gamma);
        GiNaC::ex Bc3 = B1 + B2 + B3;

        return std::make_pair(Ac3, Bc3);
    }

    std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initRangeCBF()
    {
        GiNaC::matrix d2 = x_target.sub(GiNaC::matrix{{px}, {py}});
        GiNaC::matrix R = {{GiNaC::cos(th), GiNaC::sin(th)},
                           {-GiNaC::sin(th), GiNaC::cos(th)}};
        GiNaC::matrix xt_rel = R.mul(d2);
        GiNaC::ex norm2 = GiNaC::pow(xt_rel(0,0), 2) + GiNaC::pow(xt_rel(1,0), 2);

        GiNaC::ex b4 = -norm2 + GiNaC::pow(Rs, 2);

        GiNaC::matrix grad_b4 = GiNaC::matrix(STATE_VARS, 1);
        grad_b4(0, 0) = GiNaC::diff(b4, px);
        grad_b4(1, 0) = GiNaC::diff(b4, py);
        grad_b4(2, 0) = GiNaC::diff(b4, th);
        grad_b4(3, 0) = GiNaC::diff(b4, vx);
        grad_b4(4, 0) = GiNaC::diff(b4, vy);
        grad_b4(5, 0) = GiNaC::diff(b4, w);

        GiNaC::ex lfb4 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb4 = lfb4 + grad_b4(i, 0) * f(i, 0);
        }
        
        GiNaC::matrix grad2_b4 = GiNaC::matrix(STATE_VARS, 1);
        grad2_b4(0, 0) = GiNaC::diff(lfb4, px);
        grad2_b4(1, 0) = GiNaC::diff(lfb4, py);
        grad2_b4(2, 0) = GiNaC::diff(lfb4, th);
        grad2_b4(3, 0) = GiNaC::diff(lfb4, vx);
        grad2_b4(4, 0) = GiNaC::diff(lfb4, vy);
        grad2_b4(5, 0) = GiNaC::diff(lfb4, w);
        GiNaC::ex lf2b4 = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lf2b4 = lf2b4 + grad2_b4(i, 0) * f(i, 0);
        }

        GiNaC::matrix grad_bc = GiNaC::matrix(STATE_VARS, 1);
        GiNaC::ex alpha_b4 = alpha(b4, gamma);
        grad_bc(0, 0) = GiNaC::diff(alpha_b4, px);
        grad_bc(1, 0) = GiNaC::diff(alpha_b4, py); 
        grad_bc(2, 0) = GiNaC::diff(alpha_b4, th); 
        grad_bc(3, 0) = GiNaC::diff(alpha_b4, vx); 
        grad_bc(4, 0) = GiNaC::diff(alpha_b4, vy); 
        grad_bc(5, 0) = GiNaC::diff(alpha_b4, w);
        GiNaC::ex lfb_c = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfb_c = lfb_c + grad_bc(i, 0) * f(i, 0);
        }

        GiNaC::matrix Ac4 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac4j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac4j += grad2_b4(i, 0) * g(i, j);
            }
            Ac4(0, j) = Ac4j;
        }


        GiNaC::ex B1 = lf2b4;
        GiNaC::ex B2 = lfb_c;
        GiNaC::ex B3 = alpha(lfb4 + alpha_b4, gamma);
        GiNaC::ex Bc4 = B1 + B2 + B3;

        return std::make_pair(Ac4, Bc4);

    }

    std::pair<GiNaC::matrix, GiNaC::ex> FovCBF::initVelCBF(GiNaC::ex bv)
    {
        GiNaC::matrix grad_bv = GiNaC::matrix(STATE_VARS, 1);
        grad_bv(0, 0) = GiNaC::diff(bv, px);
        grad_bv(1, 0) = GiNaC::diff(bv, py);
        grad_bv(2, 0) = GiNaC::diff(bv, th);
        grad_bv(3, 0) = GiNaC::diff(bv, vx);
        grad_bv(4, 0) = GiNaC::diff(bv, vy);
        grad_bv(5, 0) = GiNaC::diff(bv, w);
        GiNaC::ex lfbv = 0.0;
        for (int i = 0; i < STATE_VARS; i++)
        {
            lfbv = lfbv + grad_bv(i, 0) * f(i, 0);
        }
        
        GiNaC::matrix Ac_v1 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac_v1j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac_v1j += grad_bv(i, 0) * g(i, j);
            }
            Ac_v1(0, j) = Ac_v1j;
        }

        GiNaC::ex Bc_v1 = lfbv + alpha(bv, 1.0);
        return std::make_pair(Ac_v1, Bc_v1);
    }


    GiNaC::ex FovCBF::matrixSubs(GiNaC::matrix a, Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex tmp = GiNaC::subs(a, px==state(0));
        tmp = GiNaC::subs(tmp, py==state(1));
        tmp = GiNaC::subs(tmp, th==state(2));
        tmp = GiNaC::subs(tmp, vx==state(3));
        tmp = GiNaC::subs(tmp, vy==state(4));
        tmp = GiNaC::subs(tmp, w==state(5));
        tmp = GiNaC::subs(tmp, xt==target_state(0));
        tmp = GiNaC::subs(tmp, yt==target_state(1));
        return tmp;
    }

    GiNaC::ex FovCBF::valueSubs(GiNaC::ex a, Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex tmp = GiNaC::subs(a, px==state(0));
        tmp = GiNaC::subs(tmp, py==state(1));
        tmp = GiNaC::subs(tmp, th==state(2));
        tmp = GiNaC::subs(tmp, vx==state(3));
        tmp = GiNaC::subs(tmp, vy==state(4));
        tmp = GiNaC::subs(tmp, w==state(5));
        tmp = GiNaC::subs(tmp, xt==target_state(0));
        tmp = GiNaC::subs(tmp, yt==target_state(1));
        return tmp;
    }



    Eigen::VectorXd FovCBF::getSafetyConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex matrix_expr = matrixSubs(Ac_safe, state, target_state);
        Eigen::VectorXd Ac;
        Ac.resize(CONTROL_VARS);
        Ac.setZero();
        for (int i = 0; i < CONTROL_VARS; i++)
        {
            GiNaC::ex val = matrix_expr[i].evalf();
            Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        return Ac;
    }


    Eigen::VectorXd FovCBF::getLBConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex matrix_expr = matrixSubs(Ac_lb, state, target_state);
        Eigen::VectorXd Ac;
        Ac.resize(CONTROL_VARS);
        Ac.setZero();
        for (int i = 0; i < CONTROL_VARS; i++)
        {
            GiNaC::ex val = matrix_expr[i].evalf();
            Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        return Ac;
    }

    Eigen::VectorXd FovCBF::getRBConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex matrix_expr = matrixSubs(Ac_rb, state, target_state);
        Eigen::VectorXd Ac;
        Ac.resize(CONTROL_VARS);
        Ac.setZero();
        for (int i = 0; i < CONTROL_VARS; i++)
        {
            GiNaC::ex val = matrix_expr[i].evalf();
            Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        return Ac;
    }

    Eigen::VectorXd FovCBF::getRangeConstraints(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex matrix_expr = matrixSubs(Ac_range, state, target_state);
        Eigen::VectorXd Ac;
        Ac.resize(CONTROL_VARS);
        Ac.setZero();
        for (int i = 0; i < CONTROL_VARS; i++)
        {
            GiNaC::ex val = matrix_expr[i].evalf();
            Ac(i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        return Ac;
    }

    Eigen::MatrixXd FovCBF::getMaxVelContraints(Eigen::VectorXd state)
    {
        Eigen::Vector2d dummy_target;
        dummy_target.setZero();
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = matrixSubs(Ac_v1_max, state, dummy_target);
        expressions[1] = matrixSubs(Ac_v2_max, state, dummy_target);
        expressions[2] = matrixSubs(Ac_v3_max, state, dummy_target);

        Eigen::MatrixXd Acs;
        Acs.resize(CONTROL_VARS, CONTROL_VARS);
        Acs.setZero();
        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            GiNaC::ex val = expressions[i].evalf();
            Acs(i,i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        return Acs;
    }

    Eigen::MatrixXd FovCBF::getMinVelContraints(Eigen::VectorXd state)
    {
        Eigen::Vector2d dummy_target;
        dummy_target.setZero();
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = matrixSubs(Ac_v1_min, state, dummy_target);
        expressions[1] = matrixSubs(Ac_v2_min, state, dummy_target);
        expressions[2] = matrixSubs(Ac_v3_min, state, dummy_target);

        Eigen::MatrixXd Acs;
        Acs.resize(CONTROL_VARS, CONTROL_VARS);
        Acs.setZero();
        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            GiNaC::ex val = expressions[i].evalf();
            Acs(i,i) = GiNaC::ex_to<GiNaC::numeric>(val).to_double();
        }

        return Acs;
    }


    double FovCBF::getSafetyBound(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex expr = valueSubs(Bc_safe, state, target_state);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();
    
        return Bc;
    }

    double FovCBF::getRangeBound(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex expr = valueSubs(Bc_range, state, target_state);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();
    
        return Bc;
    }

    double FovCBF::getLBBound(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex expr = valueSubs(Bc_lb, state, target_state);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();
    
        return Bc;
    }

    double FovCBF::getRBBound(Eigen::VectorXd state, Eigen::VectorXd target_state)
    {
        GiNaC::ex expr = valueSubs(Bc_rb, state, target_state);
        double Bc = GiNaC::ex_to<GiNaC::numeric>(expr).to_double();
    
        return Bc;
    }

    Eigen::VectorXd FovCBF::getMaxVelBounds(Eigen::VectorXd state)
    {
        Eigen::Vector2d dummy_target;
        dummy_target.setZero();
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = valueSubs(Bc_v1_max, state, dummy_target);
        expressions[1] = valueSubs(Bc_v2_max, state, dummy_target);
        expressions[2] = valueSubs(Bc_v3_max, state, dummy_target);

        Eigen::VectorXd Bs;
        Bs.resize(CONTROL_VARS);

        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            Bs(i) = GiNaC::ex_to<GiNaC::numeric>(expressions[i]).to_double();
        }

        return Bs;
    }

    Eigen::VectorXd FovCBF::getMinVelBounds(Eigen::VectorXd state)
    {
        Eigen::Vector2d dummy_target;
        dummy_target.setZero();
        std::vector<GiNaC::ex> expressions(CONTROL_VARS);
        expressions[0] = valueSubs(Bc_v1_min, state, dummy_target);
        expressions[1] = valueSubs(Bc_v2_min, state, dummy_target);
        expressions[2] = valueSubs(Bc_v3_min, state, dummy_target);

        Eigen::VectorXd Bs;
        Bs.resize(CONTROL_VARS);

        for (int i = 0; i < CONTROL_VARS; ++i)
        {
            Bs(i) = GiNaC::ex_to<GiNaC::numeric>(expressions[i]).to_double();
        }

        return Bs;
    }

    void FovCBF::setAlpha(std::function<GiNaC::ex(GiNaC::ex, double)> newAlpha)
    {
        alpha = newAlpha;
    }
    
    



}
