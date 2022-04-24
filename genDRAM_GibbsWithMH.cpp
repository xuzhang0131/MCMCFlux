//
//  genDRAM_GibbsWithMH.cpp
//
//
//  Created by X.Z on 10/4/18.
//

#include "genDRAM_GibbsWithMH.h"
#include <iomanip>
// 'startp', give some vector as starting point,including sigma^2 and constat parameters;
// 'muprior', 'sigmaprior', assume \theta_i ~ MVN(mu_i, sigma_i^2) and mutually indep;
//
// 'indexFit', indicator for this parameter to be fitted or not with 1 or 0
//
// only have observations at 24h, but with multiple independent replicates for this time
// point
//
// mmYobs: a matrix, one row for one observation at 24h
//
//

extern const int cov_update_ind (100);
extern const int adaptint (50);
extern const double qcovadj (1e-6);

Eigen::MatrixXd ps_hist_total;

Eigen::MatrixXd chaincov;//= Eigen::MatrixXd::Zero(104,104);

void genDRAM_GibbsWithMH(int randseed, int diffterm)
{
    int i,j; // for loop
    extern const int cov_update_ind;
    extern const int adaptint;
    extern const double qcovadj;
    
    double lb_old [] = {0,0.1,0,0.1,0,0,0,0,0,50,
        0,0,0,0,0,0,0,0,0,100,
        0.5,10,0,0,0.1,0,0.2e6,100,0.1,0.2,
        0.02,0.2,1e4,10,2e4,2e4,0,0,0,0,
        1,0.05,0.1,0.01,0.5,0.1,5,1,50,1,
        1,0,
        0,0.1,0,0.1,0,0,0,0,0,50,
        0,0,0,0,0,0,0,0,0,100,
        0.5,10,0,0,0.1,0,0.2e6,100,0.1,0.2,
        0.02,0.2,1e4,10,2e4,2e4,0,0,0,0,
        1,0.05,0.1,0.01,0.5,0.1,5,1,50,1,
        1,0};
    double ub_old [] = {0.01,1,0.1,1,0.01,5e-4,1,2e-4,1e-1,500,
        10,0.01,0.1,0.001,0.2,0.02,0.03,0.01,0.05,1000,
        20,100,2,3,10,3,1e7,1000,4,10,
        1,10,1e5,100,2e5,2e5,0.25,0.01,0.5,0.3,
        50,5,8,1,10,10,100,50,500,500,
        100,3,
        0.01,1,0.1,1,0.01,5e-4,1,2e-4,1e-1,500,
        10,0.01,0.1,0.001,0.2,0.02,0.03,0.01,0.05,1000,
        20,100,2,3,10,3,1e7,1000,4,10,
        1,10,1e5,100,2e5,2e5,0.25,0.01,0.5,0.3,
        50,5,8,1,10,10,100,50,500,500,
        100,3};
    
    double p_curr [104], err_curr[n_err2];
    for(i=0; i<n_err2; ++i) err_curr[i] = startp[i+104] ;
    double lb[104];
    double ub[104];
    for(i = 0; i < 104; ++i)
    {
        p_curr[i] = startp[i];  // exlude the non-fitted parameters
        lb[i] = lb_old[i];
        ub[i] = ub_old[i];
    }
    
    //    cout << "Check0" << endl;
    double delta_com = p_curr[diffterm] - p_curr[diffterm+52];
    double lpostC = 0.0;
    double lpostN = 0.0;
    lpostLik(lpostC, p_curr,err_curr,1,diffterm);//cancer;
    lpostLik(lpostN, p_curr,err_curr,2,diffterm);//noncaner;

    
    double mybeta_err[n_err2];
    for(i=0; i<n_err2; ++i) mybeta_err[i] = beta_err[i];
    double L_currC = -2*lpostC;
    double L_currN = -2*lpostN;
    double L_curr = L_currC+L_currN;
    double L_trial;
    double L_trial2;
    
    double sp_trial;
    double qmu_c_curr;
    double qcov_c_curr;
    double tp_trial [104];
    bool bound_flag;
    double qcov2;
    double new_sp_trial;
    
    const double min_accept = 0.4;  // 0.4    0.7
    const double max_accept = 0.7;
    //extern double sd_tune;
    //    sd_tune = 2.38*2.38/104;
    
    ps_hist_total = Eigen::MatrixXd::Zero(nruns*nthinning+nburnin,104+n_err2);
    const int Cmax = 1e8;
    const double Cmin = 1e-8;
    const double Cmod = 1.1;
    //Cfactor = 0.02;
    double Cfactor [104];
    for(i=0; i<104; ++i) Cfactor[i] = 2.38*2.38/104/ 104;
    
    //    cout << "Check1" << endl;
    
    // mcmc;
    const int naccepts = 30;
    int accepts [104][naccepts];
    //accepts={0};//nan(sum(ifFit),naccepts);  //accept or not, only count for the most recent 30 ones
    int i_accepts [104];
    for(i=0; i<104; ++i) i_accepts[i] = 1; //= repmat(1,sum(ifFit),1);
    
    bool dr_flag = true;
    const int drscale = 16;
    int save_q;
    double save_r;
    const int file_save = 500;// for result data saving;
    
    //    std::default_random_engine generator;
    std::mt19937 generator;
    generator.seed(randseed*11);
    std::normal_distribution<double> distribution(0.0,1.0);
    std::uniform_real_distribution<double> distribution2(0.0,1.0);
    double randa;
    double mymin;
    double a12 = 0;
    bool qa = false;
    double a32, l2, q1, a13;
    
    
    cout << "MCMC sampling..." << endl;
    
    //int jthin = 1;
    int jcount = 1;
    int jrungo = -nburnin+1;
    
    // covariance update uses these to store previous values
    Eigen::MatrixXd CAdapt (104,104);
    double cov11;
    Eigen::MatrixXd cov12, cov22, cov_tmp;
    chaincov = Eigen::MatrixXd::Zero(104,104);
    
    for(int jruns = 0; jruns<(nruns*nthinning)+nburnin; ++jruns)
    {
        //        cout << "j" << jruns << endl;
        
        methodGibbs(jruns); //!!
        
        double accept_rate[104]; //acceptance monitor;
        
        int jGibb = 0;
        for( jGibb = 0;jGibb<52;++jGibb){
        qmu_c_curr = p_curr[jGibb];
        if(jruns<adaptint*nthinning+nburnin + cov_update_ind) qcov_c_curr = Cfactor[jGibb];
        else{
            if(ifFit.size() != 1)
            {
                reloCov(jGibb, CAdapt);//!!
                cov11 = CAdapt(CAdapt.rows()-1,CAdapt.rows()-1);
                cov12 = CAdapt.bottomLeftCorner(1,CAdapt.rows()-1);
                cov22 = CAdapt.topLeftCorner(CAdapt.rows()-1,CAdapt.rows()-1);//(1:(end-1),1:(end-1));
                cov_tmp = cov12 * (cov22.colPivHouseholderQr().solve(cov12.transpose()));
                qcov_c_curr = sd_tune[jGibb]*(cov11 - cov_tmp(0,0));
            }
            else if(ifFit.size() == 1) qcov_c_curr = sd_tune[jGibb]*chaincov(0,0);
        }
        
        sp_trial = qmu_c_curr + distribution(generator)*pow(qcov_c_curr,0.5);
        
        if(exp(sp_trial)<lb[jGibb]) sp_trial = log(exp(sp_trial) + 2*(lb[jGibb] - exp(sp_trial)));
        else if(exp(sp_trial)>ub[jGibb]) sp_trial = log(exp(sp_trial) + 2*(ub[jGibb] - exp(sp_trial)));
        
        a12 = 0;
        qa = false;
        
        bound_flag = (exp(sp_trial)>=lb[jGibb])? true:false;
        bound_flag = (exp(sp_trial)<=ub[jGibb])? true:false;
        
        for(i=0; i<104; ++i) tp_trial[i] = p_curr[i];
        
        if(bound_flag)
        {
            tp_trial[jGibb] = sp_trial;
            
            if(jGibb != diffterm){
            lpostLik(lpostC, tp_trial, err_curr,1,diffterm);
            L_trial = -2*lpostC;
            
            a12 = exp(-0.5*(L_trial - L_currC));
            randa = distribution2(generator);
            mymin = (a12<=1)?a12:1;
            qa = randa <= mymin; //% accept?
            }
            else if(jGibb == diffterm){
                tp_trial[jGibb+52] = sp_trial - delta_com;
                lpostLik(lpostC, tp_trial, err_curr,1,diffterm);
                lpostLik(lpostN, tp_trial, err_curr,2,diffterm);
                L_trial = -2*lpostC-2*lpostN;
                
                a12 = exp(-0.5*(L_trial - L_curr));
                randa = distribution2(generator);
                mymin = (a12<=1)?a12:1;
                qa = randa <= mymin;
            }
        }
        
        // update
        if(qa)   // if accept this trial
        {
            if(jGibb != diffterm){
            for(i=0; i<n_err; ++i) mybeta_err[i] = beta_err[i];
            p_curr[jGibb] = sp_trial;
            L_currC = L_trial;
                L_curr = L_currC + L_currN;
            }
            else if (jGibb == diffterm){
                for(i=0; i<n_err2; ++i) mybeta_err[i] = beta_err[i];
                p_curr[jGibb] = sp_trial;
                p_curr[jGibb+52] = sp_trial - delta_com;
                L_currC = -2*lpostC;
                L_currN = -2*lpostN;
                L_curr = L_trial;
            }
        }
        else if(!qa && dr_flag && bound_flag && jrungo>1000) //do delayed rejection
        {
            qcov2 = qcov_c_curr/drscale;
            new_sp_trial = qmu_c_curr + distribution(generator)*pow(qcov2,0.5);
            
            if(exp(new_sp_trial)<lb[jGibb]) new_sp_trial = log(exp(new_sp_trial) + 2*(lb[jGibb] - exp(new_sp_trial)));
            else if(exp(new_sp_trial)>ub[jGibb]) new_sp_trial = log(exp(new_sp_trial) + 2*(ub[jGibb] - exp(new_sp_trial)));

            bound_flag = (exp(new_sp_trial)>=lb[jGibb])? true:false;
            bound_flag = (exp(new_sp_trial)<=ub[jGibb])? true:false;
            
            if(bound_flag)
            {
                tp_trial[jGibb] = new_sp_trial;
                
                if(jGibb != diffterm){
                lpostLik(lpostC, tp_trial, err_curr,1,diffterm);
                L_trial2 = -2*lpostC;
                a32 = (exp(-0.5*(L_trial - L_trial2))<=1)? exp(-0.5*(L_trial - L_trial2)):1;
                l2 = exp(-0.5*(L_trial2-L_currC));
                q1 = exp(-0.5*(pow(sp_trial - new_sp_trial,2)/qcov_c_curr + pow(sp_trial - qmu_c_curr,2)/qcov_c_curr));
                a13 = l2*q1*(1-a32)/(1-((a12<=1)?a12:1));
                randa = distribution2(generator);
                mymin = (a13<=1)?a13:1;
                qa = randa <= mymin; //% accept?
                }
                else if(jGibb == diffterm){
                    tp_trial[jGibb+52] = new_sp_trial - delta_com;
                    lpostLik(lpostC, tp_trial, err_curr,1,diffterm);
                    lpostLik(lpostN, tp_trial, err_curr,2,diffterm);
                    L_trial2 = -2*lpostC-2*lpostN;
                    a32 = (exp(-0.5*(L_trial - L_trial2))<=1)? exp(-0.5*(L_trial - L_trial2)):1;
                    l2 = exp(-0.5*(L_trial2-L_curr));
                    q1 = exp(-0.5*(pow(sp_trial - new_sp_trial,2)/qcov_c_curr + pow(sp_trial - qmu_c_curr,2)/qcov_c_curr));
                    a13 = l2*q1*(1-a32)/(1-((a12<=1)?a12:1));
                    randa = distribution2(generator);
                    mymin = (a13<=1)?a13:1;
                    qa = randa <= mymin; //% accept?
                }
            }
            if(qa)   // if accept this trial
            {
                if(jGibb != diffterm){
                    for(i=0; i<n_err; ++i) mybeta_err[i] = beta_err[i];
                    p_curr[jGibb] = new_sp_trial;
                    L_currC = L_trial2;
                    L_curr = L_currC + L_currN;
                }
                else if (jGibb == diffterm){
                    for(i=0; i<n_err2; ++i) mybeta_err[i] = beta_err[i];
                    p_curr[jGibb] = new_sp_trial;
                    p_curr[jGibb+52] = new_sp_trial - delta_com;
                    L_currC = -2*lpostC;
                    L_currN = -2*lpostN;
                    L_curr = L_trial2;
                }
            }
        }
        
        accepts[jGibb][i_accepts[jGibb]-1] = qa? 1:0;
        accept_rate[jGibb] = std::accumulate(&accepts[jGibb][0],&accepts[jGibb][naccepts-1],0.0)/(naccepts-1);
        i_accepts[jGibb] = i_accepts[jGibb] + 1;
        if(i_accepts[jGibb]>naccepts) i_accepts[jGibb] = 1;    //only look at the accept rate of previous 30 ones
        
        // adjust scaling during burn-in
        if(jrungo<=0){
            if(accept_rate[jGibb] > max_accept && Cfactor[jGibb]*Cmod<Cmax)
                Cfactor[jGibb] = Cfactor[jGibb]*Cmod;
            else if(accept_rate[jGibb] < min_accept && Cfactor[jGibb]/Cmod>Cmin)
                Cfactor[jGibb] = Cfactor[jGibb]/Cmod;
        }
        else{ // adjust scaling after burn-in
            if(accept_rate[jGibb] > max_accept && sd_tune[jGibb]*Cmod<Cmax)
                sd_tune[jGibb] = sd_tune[jGibb]*Cmod;
            else if(accept_rate[jGibb] < min_accept && sd_tune[jGibb]/Cmod>Cmin)
                sd_tune[jGibb] = sd_tune[jGibb]/Cmod;
        }
        }// end for jGibb cancer;
        
//        L_curr  = L_currC + L_currN;
        
        //////////////////////////////////////////////////////////////
        
        
        //noncancer
        // jGibb = 52;
        for(int jGibb = 52; jGibb<104; ++jGibb)
        {
            if(ifFit[jGibb] == (52+diffterm)) {qmu_c_curr = delta_com;}
            else qmu_c_curr = p_curr[jGibb];
            
            if(jruns<adaptint*nthinning+nburnin + cov_update_ind) qcov_c_curr = Cfactor[jGibb];
            else{
                if(ifFit.size() != 1)
                {
                    reloCov(jGibb, CAdapt);//!!
                    cov11 = CAdapt(CAdapt.rows()-1,CAdapt.rows()-1);
                    cov12 = CAdapt.bottomLeftCorner(1,CAdapt.rows()-1);
                    cov22 = CAdapt.topLeftCorner(CAdapt.rows()-1,CAdapt.rows()-1);//(1:(end-1),1:(end-1));
                    cov_tmp = cov12 * (cov22.colPivHouseholderQr().solve(cov12.transpose()));
                    qcov_c_curr = sd_tune[jGibb]*(cov11 - cov_tmp(0,0));
                }
                else if(ifFit.size() == 1) qcov_c_curr = sd_tune[jGibb]*chaincov(0,0);
            }
            
            if(ifFit[jGibb] == (52+diffterm)){ //5+52
                sp_trial = p_curr[jGibb-52] - qmu_c_curr - distribution(generator)*pow(qcov_c_curr,0.5);//b12 = b11 - delta12;
            }else{sp_trial = qmu_c_curr + distribution(generator)*pow(qcov_c_curr,0.5);}
            
            
            if(exp(sp_trial)<lb[jGibb]) sp_trial = log(exp(sp_trial) + 2*(lb[jGibb] - exp(sp_trial)));
            else if(exp(sp_trial)>ub[jGibb]) sp_trial = log(exp(sp_trial) + 2*(ub[jGibb] - exp(sp_trial)));
            
            a12 = 0;
            qa = false;
            
            bound_flag = (exp(sp_trial)>=lb[jGibb])? true:false;
            bound_flag = (exp(sp_trial)<=ub[jGibb])? true:false;
            
            for(i=0; i<104; ++i) tp_trial[i] = p_curr[i];
            
            if(bound_flag)
            {
                tp_trial[jGibb] = sp_trial;
                lpostLik(lpostN, tp_trial, err_curr,2,diffterm);  //change to the partial lpostLikNoncancer();
                L_trial = -2*lpostN;
                
                a12 = exp(-0.5*(L_trial - L_currN));
                randa = distribution2(generator);
                mymin = (a12<=1)?a12:1;
                qa = randa <= mymin; //% accept?
            }
            
            // update
            if(qa)   // if accept this trial
            {
                for(i=n_err; i<n_err2; ++i) mybeta_err[i] = beta_err[i];
                p_curr[jGibb] = sp_trial;
                if(ifFit[jGibb] == (52+diffterm)) delta_com = p_curr[diffterm] - sp_trial;
                L_currN = L_trial;
            }
            else if(!qa && dr_flag && bound_flag && jrungo>1000) //do delayed rejection
            {
                qcov2 = qcov_c_curr/drscale;
                
                if(ifFit[jGibb] == (52+diffterm)){ //5+52
                    new_sp_trial = p_curr[jGibb-52] - qmu_c_curr - distribution(generator)*pow(qcov2,0.5);//b12 = b11 - delta12;
                }else{new_sp_trial = qmu_c_curr + distribution(generator)*pow(qcov2,0.5);}
                
                if(exp(new_sp_trial)<lb[jGibb]) new_sp_trial = log(exp(new_sp_trial) + 2*(lb[jGibb] - exp(new_sp_trial)));
                else if(exp(new_sp_trial)>ub[jGibb]) new_sp_trial = log(exp(new_sp_trial) + 2*(ub[jGibb] - exp(new_sp_trial)));
                
                bound_flag = (exp(new_sp_trial)>=lb[jGibb])? true:false;
                bound_flag = (exp(new_sp_trial)<=ub[jGibb])? true:false;
                
                if(bound_flag)
                {
                    tp_trial[jGibb] = new_sp_trial;
                    
                        lpostLik(lpostN, tp_trial, err_curr,2,diffterm);
                        L_trial2 = -2*lpostN;
                        a32 = (exp(-0.5*(L_trial - L_trial2))<=1)? exp(-0.5*(L_trial - L_trial2)):1;
                        l2 = exp(-0.5*(L_trial2-L_currN));
                        q1 = exp(-0.5*(pow(sp_trial - new_sp_trial,2)/qcov_c_curr + pow(sp_trial - qmu_c_curr,2)/qcov_c_curr));
                        a13 = l2*q1*(1-a32)/(1-((a12<=1)?a12:1));
                        randa = distribution2(generator);
                        mymin = (a13<=1)?a13:1;
                        qa = randa <= mymin; //% accept?
                }
                if(qa)   // if accept this trial
                {
                    for(i=n_err; i<n_err2; ++i) mybeta_err[i] = beta_err[i];
                    p_curr[jGibb] = new_sp_trial;
                    if(ifFit[jGibb] == (52+diffterm)) delta_com = p_curr[diffterm] - new_sp_trial;
                    L_currN = L_trial2;
                }
            }
            
            accepts[jGibb][i_accepts[jGibb]-1] = qa? 1:0;
            accept_rate[jGibb] = std::accumulate(&accepts[jGibb][0],&accepts[jGibb][naccepts-1],0.0)/(naccepts-1);
            i_accepts[jGibb] = i_accepts[jGibb] + 1;
            if(i_accepts[jGibb]>naccepts) i_accepts[jGibb] = 1;    //only look at the accept rate of previous 30 ones
            
            // adjust scaling during burn-in
            if(jrungo<=0){
                if(accept_rate[jGibb] > max_accept && Cfactor[jGibb]*Cmod<Cmax)
                    Cfactor[jGibb] = Cfactor[jGibb]*Cmod;
                else if(accept_rate[jGibb] < min_accept && Cfactor[jGibb]/Cmod>Cmin)
                    Cfactor[jGibb] = Cfactor[jGibb]/Cmod;
            }
            else{ // adjust scaling after burn-in
                if(accept_rate[jGibb] > max_accept && sd_tune[jGibb]*Cmod<Cmax)
                    sd_tune[jGibb] = sd_tune[jGibb]*Cmod;
                else if(accept_rate[jGibb] < min_accept && sd_tune[jGibb]/Cmod>Cmin)
                    sd_tune[jGibb] = sd_tune[jGibb]/Cmod;
            }
            
        }// end for jGibb in the noncancer group;
    
        
        for(i=0; i<104; ++i)
        {
            ps_hist_total(jruns,i) = p_curr[i];
        }
        ps_hist_total(jruns,52+diffterm) = delta_com;
        
        //error term update
        for(i=0; i<n_err2; ++i){
            std::gamma_distribution<double> mygamma(alpha_err[i], 1/mybeta_err[i]);
            err_curr[i] = 1/mygamma(generator);
            ps_hist_total(jruns,i+104) = err_curr[i];
        }
    
    
    
        lpostLik(lpostC, p_curr,err_curr,1,diffterm); // update the lpostC for cancer with new error variances;
        L_currC = -2*lpostC;
        lpostLik(lpostN, p_curr,err_curr,2,diffterm); // update the lpostN for noncancer with new error variances;
        L_currN = -2*lpostN;
        L_curr = L_currC + L_currN; // update the full lpost for next jruns round simulation;
    
    
        
        //save samples
        if(jrungo>0)
        {
            //jthin = 1;
            jcount = jcount + 1;
            //}
        }
        
        std::div_t divre = div(jruns+1,file_save);
        save_q = divre.quot;
        save_r = divre.rem;
        if(save_r == 0 & save_q > 0)
        {
            string file_name = "htsample_curr_" + to_string(diffterm) + "_" + std::to_string(save_q) + ".txt";
            ofstream myfile;
            myfile.open (file_name);
            for(i=0; i<104/2; ++i)
            {
                myfile << "par" << ifFit[i]+1 << '\t';
            }
            for(i=0; i<n_err; ++i)
            {
                myfile << "err" << i+1 << '\t';
            }
            myfile << '\n';
            
            for(i=0; i<file_save; ++i)
            {
                for(j=0; j<52; ++j) myfile << std::scientific << std::setprecision(8) << ps_hist_total(file_save*(save_q-1)+i,j) << '\t';// exp(ps[file_save*(save_q-1)+i][ifFit[j]]);
                for(j=104; j<104+n_err; ++j) myfile << std::scientific << std::setprecision(8) << ps_hist_total(file_save*(save_q-1)+i,j) << '\t';
                myfile << '\n';
            }
            
            myfile.close();
            
            string file_name2 = "htsample2_curr_" + to_string(diffterm) + "_" + std::to_string(save_q) + ".txt";
            ofstream myfile2;
            myfile2.open (file_name2);
            for(i=0; i<104/2; ++i)
            {
                myfile2 << "par" << ifFit[i]+1 << '\t';
            }
            for(i=0; i<n_err; ++i)
            {
                myfile2 << "err" << i+1 << '\t';
            }
            myfile2 << '\n';
            
            for(i=0; i<file_save; ++i)
            {
                for(j=52; j<104; ++j) myfile2 << std::scientific << std::setprecision(8) << ps_hist_total(file_save*(save_q-1)+i,j) << '\t';// exp(ps[file_save*(save_q-1)+i][ifFit[j]]);
                for(j=104+n_err; j<104+n_err2; ++j) myfile2 << std::scientific << std::setprecision(8) << ps_hist_total(file_save*(save_q-1)+i,j) << '\t';
                myfile2 << '\n';
            }
            
            myfile2.close();
        }
        
        jrungo = jrungo + 1;
    }
    
    cout << "MCMC sampling done. " << '\n';
    return;
}

void lpostLik(double& lpost, double *ptmp, double *err_curr, int myindex, int diffterm)
{
    //extern int observed [95];
    double inputp_tmp[N_tmp];
    
    int i,j; //for loop
    state_type Xd_tmp, Yd_tmp;;
    double sum_rand = 0;
    double Yode_tmp[n_Y], data_arr[34];
    double llikp_tmp = 0;
    
    if(myindex == 1)
    {
    //cancer
//    for(i=0; i<52; ++i) inputp_tmp[i] = startp[i];
    for(i=52; i<52+n_err; ++i) inputp_tmp[i] = err_curr[i-52];
    
    for(i = 0; i < 52; ++i)inputp_tmp[i] = ptmp[i];
    
    extern state_type Xd_flag;
    for(i=0; i<512; ++i) Xd_tmp[i] = Xd_flag[i];
    extern double odeParameters[57];
    for(i=0; i<(N_tmp-n_Y); i++) odeParameters[i] = exp(inputp_tmp[i]);
    integrate(PurineSynthesis, Xd_tmp, t0, tf, df);
        for(i = 0; i < Yd_tmp.size() ; ++i ){
            Yd_tmp[i] = Xd_tmp[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
        }
    //exclude unobserved data with 0
        double* data_obs = data_process(data_arr,Yd_tmp);
        for(i=0; i<n_Y; ++i) Yode_tmp[i]= log(data_obs[observed[i]]);// cout<<Yode_tmp[i]<<'\t';}
    
    //error term involved
    sum_rand = 0;
    for(i=0; i<n_err; ++i) {
        sum_rand = 0;
        for(j=0; j<3; ++j) sum_rand += pow(mmYobs[j][i]-Yode_tmp[i],2);
        beta_err[i] = beta_tmp[i]+sum_rand/2;
    }
    
    //likelihood
    for(i=0; i<3; ++i)
    {
        for(j=0; j<n_Y; ++j){
            llikp_tmp += - 0.5*pow(mmYobs[i][j]-Yode_tmp[j],2)/inputp_tmp[N_tmp-n_Y+j]
            - 0.5*log(err_curr[j]);//- 0.5*log(2.0 * M_PI)
        }
    }
    
    //prior
    for(i=0; i<52;++i)
    {
        llikp_tmp +=  - 0.5*pow(inputp_tmp[i]-muprior[i],2)/sigmaprior[i];
    }
    
    for(i=0; i<n_err; ++i)
    {
        boost::math::inverse_gamma_distribution<double> d(alpha_tmp[i], beta_tmp[i]);
        llikp_tmp += log(pdf(d,inputp_tmp[52+i]));
    }
    }// end for cancer group or myindex==1;
    
    else if (myindex == 2)
    {
    //noncancer
//    for(i=0; i<52; ++i) inputp_tmp[i] = startp[i+52];
    for(i=52; i<52+n_err; ++i) inputp_tmp[i] = err_curr[i-52+n_err];
    
    for(i = 52; i < 104; ++i)inputp_tmp[i-52] = ptmp[i];
    
    extern state_type Xd_flag2;
    for(i=0; i<512; ++i) Xd_tmp[i] = Xd_flag2[i];
            extern double odeParameters[57];
    for(i=0; i<(N_tmp-n_Y); i++) odeParameters[i] = exp(inputp_tmp[i]);
    integrate(PurineSynthesis, Xd_tmp, t0, tf, df);
        for(i = 0; i < Yd_tmp.size() ; ++i ){
            Yd_tmp[i] = Xd_tmp[i]/MoleculeNumberInOneNanoMole; // unit is nano mole
        }
        //exclude unobserved data with 0
        double* data_obs = data_process(data_arr,Yd_tmp);
        for(i=0; i<n_Y; ++i) Yode_tmp[i]= log(data_obs[observed[i]]);// cout<<Yode_tmp[i]<<'\t';}
    
    //error term involved
    sum_rand = 0;
    for(i=0; i<n_err; ++i) {
        sum_rand = 0;
        for(j=0; j<3; ++j) sum_rand += pow(mmYobs2[j][i]-Yode_tmp[i],2);
        beta_err[i+n_err] = beta_tmp[i+n_err]+sum_rand/2;
    }
    
    //likelihood
    for(i=0; i<3; ++i)
    {
        for(j=0; j<n_Y; ++j){
            llikp_tmp += - 0.5*pow(mmYobs2[i][j]-Yode_tmp[j],2)/inputp_tmp[N_tmp-n_Y+j]
            - 0.5*log(err_curr[j+n_err]);//- 0.5*log(2.0 * M_PI)
        }
    }
    
    //prior
    double delta12 = ptmp[diffterm]-ptmp[52+diffterm];
    llikp_tmp +=  - 0.5*pow(delta12-muprior[diffterm]+muprior2[diffterm],2)/sigmaprior[diffterm];
    for(i=52; i<104;++i)
    {
        if(i!= (52+diffterm)) llikp_tmp +=  - 0.5*pow(inputp_tmp[i-52]-muprior2[i-52],2)/sigmaprior[i-52];
    }
    
    for(i=n_err; i<n_err2; ++i)
    {
        boost::math::inverse_gamma_distribution<double> d(alpha_tmp[i], beta_tmp[i]);
        llikp_tmp += log(pdf(d,inputp_tmp[52+i-n_err]));
    }
    }// end for noncancer group or myindex==1;
    
    
    //    return lpost;
    lpost = llikp_tmp;
}


void methodGibbs(int& jruns)
{
    std::div_t divre = div(jruns,cov_update_ind);
    double cov_r = divre.rem;
    
    if((jruns == adaptint*nthinning+nburnin + cov_update_ind) || (cov_r == 0 && jruns > adaptint*nthinning+nburnin + cov_update_ind))
    {
        Eigen::MatrixXd ps_hist_tmp(jruns-nburnin,104);
        for(int i=0; i<(jruns-nburnin); ++i)
            for(int j=0; j< 104; ++j) ps_hist_tmp(i,j) = ps_hist_total(nburnin+i,j);
        
        Eigen::MatrixXd centered = ps_hist_tmp.rowwise() - ps_hist_tmp.colwise().mean();
        Eigen::MatrixXd cov_data = (centered.adjoint() * centered) / double(ps_hist_tmp.rows() - 1);
        //.adjoint() = .transpose() here;
        chaincov = cov_data + Eigen::MatrixXd::Identity(104, 104)*qcovadj;
    }
}

template <typename Derived>
void reloCov(int index, Eigen::MatrixBase<Derived>& CAdapt)
{
    //Eigen::MatrixXd reloM;
    int i,j;
    
    if(ifFit.size()>1 && index != 104 - 1)
    {
        Eigen::MatrixXd tmpM (chaincov), tmpM2 (chaincov);
        
        //change the rows;
        for(i=0; i<104; ++i)
            for(j=0; j<104; ++j)
            {
                if(i >= index && i <104-1) tmpM(i,j) = chaincov(i+1,j);
                else if(i == 104-1) tmpM(i,j) = chaincov(index,j);
            }
        tmpM2 = tmpM;
        //change the columns;
        for(j=0; j<104; ++j)
            for(i=0; i<104; ++i)
            {
                if(j >= index && j <104-1) tmpM(i,j) = tmpM2(i,j+1);
                else if(j == 104-1) tmpM(i,j) = tmpM2(i,index);
            }
        CAdapt = tmpM;
    }
    else if(ifFit.size()== 1 | index == 104 - 1) CAdapt = chaincov;
    
}



























