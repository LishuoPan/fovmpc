#include <vector>
#include <cbf/detail/cbf.h>


namespace cbf
{   
    // template <typename T>
    FovCBF::FovCBF(double fov, double safety_dist, double max_dist) : fov(fov), Ds(safety_dist), Rs(max_dist), px("px"), py("py"), th("th"), vx("vx"), vy("vy"), w("w"), xt("xt"), yt("yt")
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

        // ------- Safety constraint ------- //
        auto res = initSafetyCBF();
        Ac_safe = res.first;
        Bc_safe = res.second;
        // std::cout << "A1: " << Ac_safe << std::endl;
        // std::cout << "B1: " << Bc_safe.evalm() << std::endl;

        auto res2 = initBorder1CBF();
        Ac_lb = res2.first;
        Bc_lb = res2.second;
        // std::cout << "A2: " << Ac_lb << std::endl;
        // std::cout << "B2: " << bc_lb << std::endl;

        auto res3 = initBorder2CBF();
        Ac_rb = res3.first;
        Bc_rb = res3.second;
        // std::cout << "A3: " << Ac_rb << std::endl;
        // std::cout << "B3: " << bc_rb << std::endl;

        auto res4 = initRangeCBF();
        Ac_range = res4.first;
        Bc_range = res4.second;
        // std::cout << "A4: " << Ac_range << std::endl;
        // std::cout << "B4: " << Bc_range << std::endl;
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

        GiNaC::matrix Ac1 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac1j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac1j += grad2_b1(i, 0) * g(i, j);
            }
            Ac1(0, j) = Ac1j;
        }

//        // linear alpha function
//        GiNaC::ex B1 = lf2b1;
//        GiNaC::ex B2 = gamma * lfb1;
//        GiNaC::ex B3 = gamma * lfb1;
//        GiNaC::ex B4 = 0;
//        GiNaC::ex B5 = gamma * gamma * b1;

//        // cubic alpha function
        GiNaC::ex B1 = lf2b1;
        GiNaC::ex B2 = gamma * 3 * GiNaC::pow(b1, 2) * lfb1;
        GiNaC::ex B3 = gamma * GiNaC::pow(lfb1 + gamma*GiNaC::pow(b1, 3), 3);
        GiNaC::ex B4 = 0;
        GiNaC::ex B5 = 0;

        GiNaC::ex Bc1 = B1 + B2 + B3 + B4 + B5;
        
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

        GiNaC::matrix Ac2 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac2j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac2j += grad2_b2(i, 0) * g(i, j);
            }
            Ac2(0, j) = Ac2j;
        }

//        // linear alpha function
//        GiNaC::ex B12 = lf2b2;
//        GiNaC::ex B22 = gamma * lfb2;
//        GiNaC::ex B32 = gamma * lfb2;
//        GiNaC::ex B42 = 0;
//        GiNaC::ex B52 = gamma * gamma * b2;

        // cubic alpha function
        GiNaC::ex B12 = lf2b2;
        GiNaC::ex B22 = gamma * 3 * GiNaC::pow(b2, 2) * lfb2;
        GiNaC::ex B32 = gamma * GiNaC::pow(lfb2 + gamma*GiNaC::pow(b2, 3), 3);
        GiNaC::ex B42 = 0;
        GiNaC::ex B52 = 0;

        GiNaC::ex Bc2 = B12 + B22 + B32 + B42 + B52;
        
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

        GiNaC::matrix Ac3 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac3j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac3j += grad2_b3(i, 0) * g(i, j);
            }
            Ac3(0, j) = Ac3j;
        }

//        // linear alpha function
//        GiNaC::ex B13 = lf2b3;
//        GiNaC::ex B23 = gamma * lfb3;
//        GiNaC::ex B33 = gamma * lfb3;
//        GiNaC::ex B43 = 0;
//        GiNaC::ex B53 = gamma * gamma * b3;

        // cubic alpha function
        GiNaC::ex B13 = lf2b3;
        GiNaC::ex B23 = gamma * 3 * GiNaC::pow(b3, 2) * lfb3;
        GiNaC::ex B33 = gamma * GiNaC::pow(lfb3 + gamma*GiNaC::pow(b3, 3), 3);
        GiNaC::ex B43 = 0;
        GiNaC::ex B53 = 0;

        GiNaC::ex Bc3 = B13 + B23 + B33 + B43 + B53;

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

        GiNaC::matrix Ac4 =  GiNaC::matrix(1, CONTROL_VARS);
        for (int j = 0; j < CONTROL_VARS; j++) {
            GiNaC::ex Ac4j = 0.0;
            for (int i = 0; i < STATE_VARS; i++) {
                Ac4j += grad2_b4(i, 0) * g(i, j);
            }
            Ac4(0, j) = Ac4j;
        }

//        // linear alpha function
//        GiNaC::ex B14 = lf2b4;
//        GiNaC::ex B24 = gamma * lfb4;
//        GiNaC::ex B34 = gamma * lfb4;
//        GiNaC::ex B44 = 0;
//        GiNaC::ex B54 = gamma * gamma * b4;

        // cubic alpha function
        GiNaC::ex B14 = lf2b4;
        GiNaC::ex B24 = gamma * 3 * GiNaC::pow(b4, 2) * lfb4;
        GiNaC::ex B34 = gamma * GiNaC::pow(lfb4 + gamma*GiNaC::pow(b4, 3), 3);
        GiNaC::ex B44 = 0;
        GiNaC::ex B54 = 0;

        GiNaC::ex Bc4 = B14 + B24 + B34 + B44 + B54;

        return std::make_pair(Ac4, Bc4);

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

    

    



}