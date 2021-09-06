clear 
close all

n=7; %a convinent method of increascing resolution while maintaining 
    % resolution ratios and number of particles per cell

np = n.*n.*10.*1000; % number of particles
dt = 0.01;  % time step
nt = 100;  % number of time steps
length_domain = 2.5; %length of periodic element
height_domain = 0.25;

omega_shape_init = 0.25;
omega_sigma_2 = 0.25;
psi_partions_num = 20; %number of partitions in 1 state space direction
c_phi = 25;
c_t = 2;
u_mean = 0;
C_0 = 10;  % rate parameter for velocity change
stochiometry = [1,1];
activity = [1,1];
B=(1+1.5*C_0); %for reducing the langevin to a standard form - part of boundary conditions implementaion from Erban and Chapman 06
x_res=n*10;%number of cells in x dim
y_res=n;

%nu and Re now used for time scale, which is turned off
% nu = 1.5*10^(-5);
% Re = 46;

n_tests = 50;

%using gamma priors, as these are scale parameters - so would be expected
%to be log normal. The gamma approximates this while having shorter tails,
%so reduces required sampling
bc_rate_prior = makedist('Gamma',1/1,1);
bc_equib_prior = makedist('Gamma',1/1,1);
omega_mean_prior = makedist('Gamma',1/1,1);
turb_e_init_prior = makedist('Gamma',1.5/0.25,0.25);
coeffs = zeros(n_tests,4);
coeffs(:,1)=random(bc_rate_prior,n_tests,1);
coeffs(:,2)=random(bc_equib_prior,n_tests,1);
coeffs(:,3)=random(omega_mean_prior,n_tests,1);
coeffs(:,4)=random(turb_e_init_prior,n_tests,1);

for index_output=1:n_tests
    bc_rate=coeffs(index_output,1);
    bc_equib=coeffs(index_output,2);
    omega_mean=coeffs(index_output,3);
    turb_e_init=coeffs(index_output,4);
    T_omega = 1/(omega_mean); %approximation to frequency timescale
    
    xp = length_domain*rand(np,nt+1); % position of particle in x-direction
    yp = height_domain*rand(np,nt+1); % position of paticle in y-direction
    phip = zeros(2, np, nt+1);%scalar concentration at these points

    %triple delta start
    for i=1:np
        a = rand();
        if a<1/3
            delta_neg1(i) = true;
            delta_plus1(i) = false;
            delta_phi2(i) = false;
        elseif a>2/3
            delta_neg1(i) = false;
            delta_plus1(i) = true;
            delta_phi2(i) = false;
        else
            delta_neg1(i) = false;
            delta_plus1(i) = false;
            delta_phi2(i) = true;
        end
    end
    
    phip(1,[delta_neg1],1) = -sqrt(3)/2+0.001*randn(sum(delta_neg1),1);
    phip(2,[delta_neg1],1) = -0.5+0.001; %pdf can't find zeros
    phip(1,[delta_plus1],1) = sqrt(3)/2+0.001*randn(sum(delta_plus1),1);
    phip(2,[delta_plus1],1) = -0.5+0.001; %pdf can't find zeros
    phip(1,[delta_phi2],1) = 0.001;%pdf can't find zeros
    phip(2,[delta_phi2],1) = 1.0+0.001*randn(sum(delta_phi2),1); 

    % top/bottom rows

%     delta_phi1 =  yp(:,1)>0.5*height_domain;
%     delta_phi2 =  logical(1-delta_phi1);
%     phip(1,[delta_phi2],1) = 0.001;%pdf can't find zeros
%     phip(2,[delta_phi2],1) = 1.0+0.001*randn(sum(delta_phi2),1); 
%     phip(2,[delta_phi1],1) = 0.001;%pdf can't find zeros
%     phip(1,[delta_phi1],1) = 1+0.001*randn(sum(delta_phi1),1); %+1;

    omegap = zeros(np,nt+1);%turbulence frequency, might want to link to rw varience
    gamran = makedist('Gamma',omega_mean/omega_shape_init,omega_shape_init);
    omegap(:,1) = random(gamran,np,1);

    %boundaries for psp cells in format x_min, x_max, y_min, y_max
    x_edges = linspace(0,length_domain,x_res+1);
    y_edges = linspace(0,height_domain,y_res+1);
    func = @(i) [x_edges(floor(i/y_res)+1),x_edges(floor(i/y_res)+2),y_edges(mod(i,y_res)+1),y_edges(mod(i,y_res)+2)];
    psp_cell_boundaries = zeros(x_res*y_res,4);
    for i=0:x_res*y_res-1
        psp_cell_boundaries(i+1,:)=func(i);
    end

    %variables for mean normalisation
    phi_pos_mean = zeros(2,1);
    phi_neg_mean = zeros(2,1);  
    phi_mean = zeros(2,1);
    a = zeros(2,np);
    phi_mean_store = zeros(2,length(psp_cell_boundaries),nt);
    c_mean_diff = zeros(np,nt);

    cond_diff = zeros(2,psi_partions_num,psi_partions_num,length(psp_cell_boundaries),nt+1);

    phi_pm = zeros(2, np); %phi plus and minus indicies to find centres for phi

    f_phi = zeros(psi_partions_num,psi_partions_num,length(psp_cell_boundaries),nt+1); 

    psi_1 = linspace(-1.2, 1.2,psi_partions_num+1);%defines boundaries for cells in psi space
    psi_2 = linspace(-1.2, 1.2,psi_partions_num+1);

    count=0;
    %populating the list
    for bound_i =1:length(psp_cell_boundaries)
        x_bounds = (xp(:,1)>=psp_cell_boundaries(bound_i,1)).*(xp(:,1)<psp_cell_boundaries(bound_i,2));
        y_bounds = (yp(:,1)>=psp_cell_boundaries(bound_i,3)).*(yp(:,1)<psp_cell_boundaries(bound_i,4));
        both_bounds = x_bounds.*y_bounds;
        cell_points = find(both_bounds);
        for j=1:sum(both_bounds)
            i=cell_points(j);
            count=0;
            while true
                count = count+1;
                phi_pm_test = randsample(cell_points,2,false);
                %applying the dot product test condition
                test_dot = dot((phip(:,phi_pm_test(1),1)-phip(:,i,1)),(phip(:,phi_pm_test(2),1)-phip(:,i,1))) ;
                if test_dot<=0
                    phi_pm(:,i) = phi_pm_test;
                    break
                end
            end
        end
    end

    %ititialising f_phi/stored data
    i=0;
    for bound_i =1:length(psp_cell_boundaries)
            x_bounds = (xp(:,i+1)>=psp_cell_boundaries(bound_i,1)).*(xp(:,i+1)<psp_cell_boundaries(bound_i,2));
            y_bounds = (yp(:,i+1)>=psp_cell_boundaries(bound_i,3)).*(yp(:,i+1)<psp_cell_boundaries(bound_i,4));
            both_bounds = x_bounds.*y_bounds;
            cell_points = find(both_bounds);
            for psi1_i=1:psi_partions_num
                for psi2_i=1:psi_partions_num
                    %histogram of phis at point
                    cond_phi1 = (phip(1,cell_points,i+1) > psi_1(psi1_i)).*(phip(1,cell_points,i+1) <=  psi_1(psi1_i+1));
                    cond_phi2 = (phip(2,cell_points,i+1) > psi_2(psi2_i)).*(phip(2,cell_points,i+1) <=  psi_2(psi2_i+1));
                    cond_both_phi = logical(cond_phi1.*cond_phi2);
                    f_phi(psi1_i,psi2_i,bound_i,i+1) = sum(cond_both_phi);
                end
            end
            %normalising f_phi
            f_phi(:,:,bound_i,i+1) = f_phi(:,:,bound_i,i+1)/sum(both_bounds);
    %         if sum(f_phi(:,:,bound_i,i+1))==0 
    %             disp('failure');disp(bound_i);disp(i+1);disp(sum(both_bounds))
    %         end
    end

    %time stamp until new p/m bound found, needed to ensure particles are
    %decorreltaed
    t_decorr_p = 1/(c_t*omegap(phi_pm(1,:),1));
    t_decorr_m = 1/(c_t*omegap(phi_pm(2,:),1));

    %used for lagvagian 
    uxp = zeros(np,1);
    uyp = zeros(np,1);

    dphi_for_plot = zeros(2,nt);
    phi_c = zeros(2,np,nt);
    system_N=zeros(nt+1,1);
    system_N(1)=np;

    k=turb_e_init;
    %fultuctaing velocities, initialised to have correct energy
    uxp = randn(np,1)*sqrt(2/3*k); 
    uyp = randn(np,1)*sqrt(2/3*k);k=zeros(np,1);
    true_time=0;
    for i = 1:nt
        time_normalised = dt*i;

        k_mean=0;
        %lagvagian approach
        for bound_i =1:length(psp_cell_boundaries)
            x_bounds = (xp(:,i+1)>=psp_cell_boundaries(bound_i,1)).*(xp(:,i+1)<psp_cell_boundaries(bound_i,2));
            y_bounds = (yp(:,i+1)>=psp_cell_boundaries(bound_i,3)).*(yp(:,i+1)<psp_cell_boundaries(bound_i,4));
            both_bounds = x_bounds.*y_bounds;
            cell_points = find(both_bounds);
    %         If terms that change energy are included, this needs recalculated
    %         - if not don't update to improve stability
            k(cell_points)=0.5*(mean(uxp(cell_points).^2)+mean(uyp(cell_points).^2))*3/2;%turb_e_init; %final term is a dimesionality correction
            k_mean = k_mean+0.5*(mean(uxp(cell_points).^2)+mean(uyp(cell_points).^2))*3/2;
            %fultuctaing velocities, as implemented in Meyer 2010 or code for Meyer
            %and Jenny 2006
    %         D(cell_points)= C_0*k/(B^2);
            uxp(cell_points) = uxp(cell_points)+(-0.5*B*omega_mean*uxp(cell_points))*dt+randn(length(cell_points),1).*sqrt(C_0.*k(cell_points).*omega_mean.*dt); 
            uyp(cell_points) = uyp(cell_points)+(-0.5*B*omega_mean*uyp(cell_points))*dt+randn(length(cell_points),1).*sqrt(C_0.*k(cell_points).*omega_mean.*dt); 
        end 
        xp(:,i+1) = xp(:,i) + (u_mean+uxp(:))*dt; % random walk in x-direction
        yp(:,i+1) = yp(:,i) + uyp(:)*dt; % random walk in y-direction
        k_mean = k_mean/length(psp_cell_boundaries);

        if false %haven't been able to get this to work
            %calculating true time scale
            C_ts = 1.5;
            p0_ts = 2;
            beta_ts = 5.2;

            %apparently these are reynolds dependent - doesn't give a way to find
            %them though
            cl_ts = 3.948;
            cn_ts=0.4304;

            int_eps = 0.0001; %an required to compute integral, as integrad(0)=NaN

            n_kol = (nu^3*k_mean/omega_mean)^2;
            L=k_mean^(2.5)/omega_mean;
            length_scale_u = sqrt(2*k_mean);
            f_l = @(x) (x./sqrt(x.^2+cl_ts)).^(5./3+p0_ts);
            f_n = @(x) exp(-beta_ts./((x.^4+cn_ts.^4).^(0.25)-cn_ts));
            integrand=@(x) x.^(-8/3).*f_l(x.*L/n_kol).*f_n(x);
            int_ts = 3.*pi.*C_ts.*0.25.*(20/3).^(1.25).*Re.^(-2.5).*integral(integrand,int_eps,Inf);
            e_turnover = int_ts*L/length_scale_u;
            correction=e_turnover;
            true_time = true_time + correction*dt;
            disp(true_time)
        end 

        % Reflection particles at boundaries

        % Reflection at upper boundary y>height_domain
        %doing closed on top open on bottom, as cell detection is open on top,
        %closed on bottom
        mag= find(yp(:,i+1)>=height_domain); % index of particle with yp>height_domain
        dim_mag = size(mag); % dimension of array "mag"

        y_mag_succ = yp(mag,i+1); % yp at time t+1 corresponding to the index "mag"

        V1 = height_domain*ones(dim_mag); % unitary array with size equat to "mag"

        ypr_mag = V1*2 -y_mag_succ ; % yp at time t+1 of the reflected particle

        yp([mag],i+1)= ypr_mag; %replacement of yp>1 with yp of reflected particle
        uyp([mag]) = -uyp([mag]); %reflecting velocity

        %is absorbed at bc?
        Q=prod(stochiometry.'.*phip(:,mag),1);
        rate = -bc_rate *(1-bc_equib*Q).';
        abs_ratio = stochiometry.*rate;
        %ratio to react at edge
        P = abs_ratio.*sqrt(B.*pi./(C_0.*k([mag])));
        xi = rand(dim_mag(1),2);
        phip(1,[mag(P(:,1)>xi(:,1))],i) = 0.0001;
        phip(2,[mag(P(:,2)>xi(:,2))],i) = 0.0001;

        % Reflection at lower boudary y<0
        neg= find(yp(:,i+1)<=0); % index of particle with yp<0
        dim_neg = size(neg); % dimension of array "mag"

        y_neg_succ = yp(neg,i+1); % dimension of array "mag"

        ypr_neg = abs(y_neg_succ);  % yp at time t+1 of the reflected particle
        yp([neg],i+1) = ypr_neg; %replacement of yp<0 with yp of reflected particle
        uyp([neg]) = -uyp([neg]); %reflecting velocity
        %is absorbed at bc?
        Q=prod(stochiometry.'.*phip(:,neg),1);
        rate = -bc_rate *(1-bc_equib*Q).';
        abs_ratio = stochiometry.*rate;
        P = abs_ratio.*sqrt(B.*pi./(k([neg]).*C_0));
        xi = rand(dim_neg(1),2);
        phip(1,[neg(P(:,1)>xi(:,1))],i) = 0.0001;
        phip(2,[neg(P(:,2)>xi(:,2))],i) = 0.0001;

        %period BC at end of periodic cell
        end_indicies = find(xp(:,i+1)>=length_domain); % index of particle with xp>length

        end_x = xp(end_indicies,i+1);
        xpr_end = end_x - length_domain; %shifting particles back to begining
        xp([end_indicies],i+1) = xpr_end; %replacing x coords

        %resetting concentrations to a bc to simulate continuous release of
        %scalar at edge of domain use comments to turn on/off
        delta_x = (yp([end_indicies],i+1)>0.5*height_domain);
        phip(1,end_indicies,i) = delta_x.*(1.0+0.001*randn(length(end_indicies),1))+ ...
            (1-delta_x)*0.001;
        phip(2,end_indicies,i) = delta_x*0.001+ ...
            (1-delta_x).*(1+0.001*randn(length(end_indicies),1));

        %resetting energy of inputted particles
        uxp(end_indicies) = randn(length(end_indicies),1)*sqrt(2/3*turb_e_init); 
        uyp(end_indicies) = randn(length(end_indicies),1)*sqrt(2/3*turb_e_init);


    %can turn off start to end periodicity if mean flow is suffiecent
        start_indicies = find(xp(:,i+1)<=0); % index of particle with xp>length
        
        start_x = xp(start_indicies,i+1);
        xp([start_indicies],i+1) = -start_x; %reflecting
        uxp([start_indicies]) = -uxp([start_indicies]); 

        %handeling php
        %E-M solver for omega
        dw = sqrt(dt).*randn(np,1);%random draws
        omegap(:,i+1) = omegap(:,i)-(omegap(:,i)-omega_mean)./T_omega.*dt + sqrt(2*omega_sigma_2*omegap(:,i).*omega_mean./T_omega).*dw;
        omegap(:,i+1) = omegap(:,i+1).*(omegap(:,i+1)>0); %enforcing positivness

        %stepping the decorrelation times
        t_decorr_p = t_decorr_p-dt;
        t_decorr_m = t_decorr_m-dt;

        %split into cells, compute centres/targets, run ODE step
        w=ones(np,1);
        for bound_i =1:length(psp_cell_boundaries)
            %test if in psp cell
            x_bounds = (xp(:,i+1)>=psp_cell_boundaries(bound_i,1)).*(xp(:,i+1)<psp_cell_boundaries(bound_i,2));
            y_bounds = (yp(:,i+1)>=psp_cell_boundaries(bound_i,3)).*(yp(:,i+1)<psp_cell_boundaries(bound_i,4));
            both_bounds = logical(x_bounds.*y_bounds);
            cell_points = find(both_bounds);
            %test if re-allocation of bounds is needed
            phi_plus_in = ismember(phi_pm(1,[cell_points]),cell_points);
            phi_minus_in = ismember(phi_pm(2,[cell_points]),cell_points);
            phi_pm_in = phi_plus_in.*phi_minus_in;
            phi_pm_in_index = find(phi_pm_in);
            phi_p_diff = (phip(:,phi_pm(1,cell_points([phi_pm_in_index])),i)-phip(:,phi_pm(1,cell_points([phi_pm_in_index])),i));
            phi_m_diff = (phip(:,phi_pm(2,cell_points([phi_pm_in_index])),i)-phip(:,phi_pm(1,cell_points([phi_pm_in_index])),i));
            test_dot = dot(phi_p_diff,phi_m_diff, 1);
            both_in_succ = find(test_dot>=0);
            redo_pairing = logical(1-ismember(both_in_succ,cell_points));
            loop_points = cell_points(redo_pairing);
            loop_len = sum(redo_pairing);
            for j=1:loop_len
                count=0;
                check_index=loop_points(j);
                while true
                    count = count+1;
                    phi_pm_test = randsample(cell_points,2,false);
                    %applying the dot product test condition
                    test_dot = dot((phip(:,phi_pm_test(1),1)-phip(:,check_index,1)),(phip(:,phi_pm_test(2),1)-phip(:,check_index,1))) ;
                    if test_dot<=0
                        phi_pm(:,check_index) = phi_pm_test;    
                        t_decorr_p(check_index) = 1/(c_t*omegap(phi_pm(1,check_index),1));
                        t_decorr_m(check_index) = 1/(c_t*omegap(phi_pm(2,check_index),1));
                        if count>np
                            disp('long count')
                        end
                        break
                    end
                end
            end
            if any(t_decorr_p<0)
                redo_pos = cell_points(t_decorr_p(cell_points)<0);
                for j=1:length(redo_pos)
                    index=redo_pos(j);
                    count=0;
                    while true
                        count=count+1;
                        phi_p_test = randsample(cell_points,1,false);
                        test_dot = dot((phip(:,phi_p_test,1)-phip(:,index,1)),(phip(:,phi_pm(2,index),1)-phip(:,index,1)));
                        if test_dot<=0
                            phi_pm(1,index) = phi_p_test;
                            t_decorr_p(index) = 1/(c_t*omegap(phi_pm(1,index),1));
                            break
                        end
                        if count>np
                            disp('long count')
                        end
                    end
                end
            end
            if any(t_decorr_m<0)
                redo_neg = cell_points(t_decorr_m(cell_points)<0);
                for j=1:length(redo_neg)
                    index=redo_neg(j);
                    count=0;
                    while true
                        count=count+1;
                        phi_m_test = randsample(cell_points,1,false);
                        test_dot = dot((phip(:,phi_m_test,1)-phip(:,index,1)),(phip(:,phi_pm(1,index),1)-phip(:,index,1)));
                        if test_dot<=0
                            phi_pm(2,index) = phi_p_test;
                            t_decorr_m(index) = 1/(c_t*omegap(phi_pm(2,index),1));
                            break
                        end
                        if count>np
                            disp('long count')
                        end
                    end
                end
            end
        end

        phi_c(:,:,i) = 0.5.*(phip(:,phi_pm(1,:),i)+phip(:,phi_pm(2,:),i));
        diffusion = c_phi.*0.5.*omegap(:,i).'.*(phip(:,:,i)-phi_c(:,:,i));
        reaction = 0;%[-1;-1].*(phip(:,:,i)./(phip(1,:,i).^2+phip(2,:,i).^2));% body reaction
        dphi = (-diffusion + reaction)*dt;
    %     ensuring mean 0 change
    %     generating a random orthonormal basis
    %     is 2-d so genrate a random unit vector from an angle and proceed based
    %     on that
        angle = 2*pi*rand(1);
        e_1 = [cos(angle);sin(angle)];
        handedness = randsample([-1,1],1);%randomly choose betwen left or right handed system
        e_2 = handedness*[[0,1];[-1,0]]*e_1;
        T = [e_1,e_2];%coord transform matrix
        dphi = T\dphi;%transform to new coords
    %     adjusting the mean change
        for bound_i =1:length(psp_cell_boundaries)
            x_bounds = (xp(:,i+1)>=psp_cell_boundaries(bound_i,1)).*(xp(:,i+1)<psp_cell_boundaries(bound_i,2));
            y_bounds = (yp(:,i+1)>=psp_cell_boundaries(bound_i,3)).*(yp(:,i+1)<psp_cell_boundaries(bound_i,4));
            both_bounds = x_bounds.*y_bounds;
            cell_points = find(both_bounds);
            for phi_i=1:2
                cell_points_pos = dphi(phi_i,cell_points)>0;
                phi_mean(phi_i) = mean(dphi(phi_i,cell_points));
                phi_pos_mean(phi_i) = mean(cell_points_pos.*dphi(phi_i,cell_points));
                phi_neg_mean(phi_i) = mean((1-cell_points_pos).*dphi(phi_i,cell_points));
                if phi_mean(phi_i)>0
                    a(phi_i,cell_points)=-cell_points_pos*(phi_neg_mean(phi_i)/phi_pos_mean(phi_i)) + (1-cell_points_pos);
                else
                    a(phi_i,cell_points)=-(1-cell_points_pos)*(phi_pos_mean(phi_i)/phi_neg_mean(phi_i)) + cell_points_pos;
                end
            end
            phi_mean_store(:,bound_i,i) = phi_mean;
            temp_ = phi_c(:,cell_points)-phi_mean;
            c_mean_diff(cell_points,i) = sqrt(temp_(1,:).^2+temp_(2,:).^2);
        end
        dphi = a.*dphi;
        dphi = T*dphi;%return to old coords

        phip(:,:,i+1) = phip(:,:,i)+dphi;

        %evaluating gamma(cond_diff) and f_phi
        for bound_i =1:length(psp_cell_boundaries)
            x_bounds = (xp(:,i+1)>=psp_cell_boundaries(bound_i,1)).*(xp(:,i+1)<psp_cell_boundaries(bound_i,2));
            y_bounds = (yp(:,i+1)>=psp_cell_boundaries(bound_i,3)).*(yp(:,i+1)<psp_cell_boundaries(bound_i,4));
            both_bounds = x_bounds.*y_bounds;
            cell_points = find(both_bounds);
            for psi1_i=1:psi_partions_num
                for psi2_i=1:psi_partions_num
                    %histogram of phis at point
                    cond_phi1 = (phip(1,cell_points,i+1) >= psi_1(psi1_i)).*(phip(1,cell_points,i+1) < psi_1(psi1_i+1));
                    cond_phi2 = (phip(2,cell_points,i+1) >= psi_2(psi2_i)).*(phip(2,cell_points,i+1) < psi_2(psi2_i+1));
                    cond_both_phi = logical(cond_phi1.*cond_phi2);
                    f_phi(psi1_i,psi2_i,bound_i,i+1) = sum(cond_both_phi);
                end
            end
            %normalising f_phi
            f_phi(:,:,bound_i,i+1) = f_phi(:,:,bound_i,i+1)/sum(both_bounds);
            system_N(i+1) = system_N(i+1)+sum(both_bounds);
        end        
    end
    [mesh_psi1,mesh_psi2 ]= meshgrid(0.5*(psi_1(1:end-1)+psi_1(2:end)),0.5*(psi_2(1:end-1)+psi_2(2:end)));
    % meshgrid output the transpose of how f_phi is labeled
    mesh_psi1=mesh_psi1.';
    mesh_psi2 = mesh_psi2.';
    % This gives equal weight to all cells, even those with less particles
    spacially_ave_f = squeeze(mean(f_phi(:,:,:,end),3));
    mesh_phis_ = cat(1,reshape(mesh_psi1,1,psi_partions_num,psi_partions_num),reshape(mesh_psi2,1,psi_partions_num,psi_partions_num));
    spacially_ave_f_ = reshape(spacially_ave_f,1,psi_partions_num,psi_partions_num);
    means(:,index_output) = squeeze(sum(mesh_phis_.*spacially_ave_f_,[2,3]));
    vars(:,index_output) = sum((mesh_phis_-means(:,index_output)).*(mesh_phis_-means(:,index_output)).*spacially_ave_f_,[2,3]);
    skews(:,index_output) = sum((mesh_phis_-means(:,index_output)).^3.*spacially_ave_f_,[2,3]);
    kurts(:,index_output) = sum((mesh_phis_-means(:,index_output)).^4.*spacially_ave_f_,[2,3]);
end

%note, just refering to central moments here, not standardiesd
%standersiation is required as part of the data analysis
save('means.mat','means')
save('vars.mat','vars')
save('skews.mat','skews')
save('kurts.mat','kurts')
save('coeffs.mat','coeffs')
clear
close 'all'