function CantileverRodSMC
rng(1, 'twister' );
L =0.4 ; N = 100; E1 = 210e9;r = 0.001 ;
rho = 8000;vstar =[0;0;1];vsstar = [0;0;0];
s1 = 20;k = 20;del_s = L/N;xi = 1;
g = 10*[0;0;0];
del_t = 0.005 ;
STEPS = 300 ;
xd = [0.115;0;0.38];

alphaf = 70/0.05;mu = 0.05;
% alphaf = 15;mu = 0.15;
fadd = [0;0;0];
F_tip = [0;0;0];
M_tip = [0;0;0];
p0 = [0;0;0];
h0 = [1;0;0;0];
q0 = [0;0;0];
w0 = [0;0;0];

xd = [0.3;0;0.38];
x2d = [0;0;0];
x3d = [0;0;0];


flist = [];
A= pi*r^2 ;G1 = E1/(2*(1+0.3));ds =  L/(N-1);
J=diag([pi*r^4/4 pi*r^4/4 pi*r^4/2]);
Kse = diag([G1*A;G1*A;G1*A]);
Kse_inv = Kse^-1;
Kbt = diag([E1*J(1,1), E1*J(2,2), G1*J(3,3)]);
Kbt_inv = Kbt^-1;
Kse_vstar = Kse*vstar;
c0 = 1.5/del_t ; c1 = -2/del_t;c2 = 0.5/del_t;


rhoA = rho*A; rhoAg = rho*A*g; rhoJ = rho*J;
y = [zeros(2,N);linspace(0,L,N);zeros(16,N)];
z = [zeros(2,N); ones(1,N); zeros(3,N)];
% fdist = [0;0;0];
y_prev = y; z_prev = z;
visualize();
G = zeros(6,1);

xend = [];
xend2 = [];
xdlist = [];
fdlist = [];
elist = [];
tlist = [];
position_x_list = [];
position_y_list = [];

for i = 2 : STEPS
    tlist = [tlist (i-2)*del_t];
    fdist = [rho*A*200*(rand()-0.5);0;0];
    % fdist = [0;0;0];
    yh = c1*y + c2*y_prev; zh = c1*z + c2*z_prev;
    y_prev = y; z_prev = z;
    yh_int = 0.5*(yh(:,1:end-1) + yh(:,2:end));
    zh_int = 0.5*(zh(:,1:end-1) + zh(:,2:end));
    G = fsolve((@(G) getResidual(G)), G);
    qL = y(1:3,N);
    fadd = getcontrol();
    visualize();
    xend = [xend y(1,end)];
    xend2 = [xend2 y(3,end)];
    t = del_t*i;
    xdlist = [xdlist xd(1)];
    flist = [flist fadd(1)];
    fdlist = [fdlist fdist(1)];
    position_x_list = [position_x_list ; y(3,:)];
    position_y_list = [position_y_list ; y(1,:)];
end

%% Visulize
figure(2)
plot(linspace(0,STEPS*del_s,STEPS-1),xend);
hold on
plot(linspace(0,STEPS*del_s,STEPS-1),xdlist);
figure(3)
plot(flist*del_s/(rho*A));
figure(4);
plot(fdlist*del_s/(rho*A));
figure(5);
surf(repmat(tlist', 1, N), position_x_list, position_y_list);
timesurt = repmat(tlist', 1, N);
shading interp; 
colormap('parula'); 
caxis([min(position_y_list(:)) max(position_y_list(:))]); 
view(-30, 30); 
hold on;
sparseIndex = 1:10:N;
sparseIndex_t = 1:2:STEPS; 
sparseIndex_t = [sparseIndex_t,STEPS-1];
sparseIndex = [sparseIndex,N];
sparseTList = timesurt(:, sparseIndex);
sparseXData = position_x_list(:, sparseIndex);
sparseYData = position_y_list(:, sparseIndex);
sparseTList = sparseTList(sparseIndex_t,:);
sparseXData = sparseXData(sparseIndex_t,:);
sparseYData = sparseYData(sparseIndex_t,:);
mesh(sparseTList, sparseXData, sparseYData, 'FaceColor', 'none', 'EdgeColor', 'k');
xlabel('Time [s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('X [m]', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
zlabel('Z [m]', 'FontName', 'Times New Roman', 'FontWeight', 'bold');

    function visualize()
        figure(1)
        plot(y(3,:),y(1,:)); title ( 'Cantilever Rod');
        xlabel('z (m)'); ylabel('x (m)');
        axis([0 1.1*L -L L]);
        grid on; daspect([1 1 1]); drawnow;
    end
    function res = sat(input)
        limit = 150;
        res = zeros(3,1);
        for iter = 1:size(input)
            res(iter) = input(iter);
            if res(iter)>limit
                res(iter) = limit;
            elseif res(iter) <-limit
                res(iter)=limit;
            end
        end
    end
    function E = getResidual(G)
        n0 = G(1:3);m0 = G(4:6);
        y(:,1) = [p0;h0;n0;m0;q0;w0];
        for j = 1:N-1
            yj = y(:,j);
            yhj_int = yh_int(:,j);
            E1 = 210e9+10e9*j/N-1;
            G1 = E1/(2*(1+0.3));
            Kse = diag([G1*A;G1*A;G1*A]);
            Kse_inv = Kse^-1;
            Kbt = diag([E1*J(1,1), E1*J(2,2), G1*J(3,3)]);
            Kbt_inv = Kbt^-1;
            Kse_vstar = Kse*vstar;
            [k1,z(:,j)]=ODE(yj,yh(:,j),zh(:,j),j);
            [k2,~]=ODE(yj + k1*ds/2,yhj_int, zh_int(:,j),j);
            [k3,~]=ODE(yj + k2*ds/2,yhj_int, zh_int(:,j),j);
            [k4,~]=ODE(yj + k3*ds,yh(:,j + 1), zh(:,j + 1),j);
            y(:,j + 1) = yj + ds*(k1 + 2*(k2 + k3) + k4)/6;
        end
        nL = y(8:10,N); mL = y(11:13,N);
        E = [F_tip-nL; M_tip-mL];
    end
    function [ys,z] = ODE(y,yh,zh,j)
        h = y(4:7); n = y(8:10); m = y(11:13);
        q = y(14:16); w = y(17:19);
        vh = zh(1:3); uh = zh(4:6);
        h1=h(1); h2=h(2); h3=h(3); h4=h(4);
        R = eye(3) + 2/(h'*h) *...
            [-h3^2-h4^2, h2*h3-h4*h1, h2*h4 + h3*h1;  h2*h3 + h4*h1, -h2^2-h4^2, h3*h4-h2*h1; h2*h4-h3*h1, h3*h4 + h2*h1, -h2^2-h3^2];
        v=Kse_inv*(R'*n + Kse_vstar);
        u=Kbt_inv*(R'* m);
        z=[v;u];
        yt = c0*y + yh; zt = c0*z + zh;
        vt = zt(1:3);
        ut = zt(4:6);
        qt = yt(14:16);
        wt = yt(17:19);
        f = rhoAg+fadd+fdist;
        ps = R*v;
        ns = rhoA*R*(cross(w,q) + qt) - f;
        ms = R*(cross(w,rhoJ*w) + rhoJ*wt)-cross(ps,n);
        qs = vt - cross(u,q) + cross(w,v);
        ws = ut - cross(u,w);
        hs = [ 0, -u(1), -u(2), -u(3);
            u(1), 0, u(3), -u(2);
            u(2), -u(3), 0, u(1);
            u(3), u(2), -u(1), 0 ] * h/2;
        ys = [ps;hs;ns;ms;qs;ws];
    end
    function f = getcontrol()
        f = zeros(3,1);
        % control parameter
        Ec = 210e9;
        Gc = Ec/(2*(1+0.3));
        Ksec = diag([Gc*A;Gc*A;Gc*A]);
        f_c = 80*rhoA;% with model error
        pL = y(1:3,end);hL = y(4:7,end);qL = y(14:16,end);
        wL = y(17:19,end);uL = z(4:6,end);vL = z(1:3,end);
        vLh1 = z(1:3,end-1);vLh2 = z(1:3,end-2);
        h1=hL(1); h2=hL(2); h3=hL(3); h4=hL(4);
        R = eye(3) + 2/(hL'*hL) *...
            [-h3^2-h4^2, h2*h3-h4*h1, h2*h4 + h3*h1;  h2*h3 + h4*h1, -h2^2-h4^2, h3*h4-h2*h1; h2*h4-h3*h1, h3*h4 + h2*h1, -h2^2-h3^2];
        x1 = pL;
        x2 = R*qL;
        vsL = 1.5*vL/del_s-2*vLh1/del_s+0.5*vLh2/del_s;
        ksev = Ksec*(vL-vstar);
        beta = Ksec*(vsL-vsstar)+cross(uL,ksev)-rho*A*cross(wL,qL);
        Omega = R*(cross(wL,qL))+R/(rho*A)*(beta);

        e1 = x1-xd;
        e1 = e1(1);
        e2 = x2-x2d;
        e2 = e2(1);

        elist = [elist e1];
        if length(elist)>3
            s = s1*e1(1)+e2(1)+alphaf*mu*glfdiff(elist,tlist,mu);
            f(1) = rho*A*(x3d(1)-s1*(x2(1)-x2d(1))-xi*sign(s)-k*s-Omega(1)-alphaf*mu*glfdiff(elist,tlist,mu-1))-f_c*sign(s)-rhoA*g(1);

        else
            s = s1*e1(1)+e2(1);
            f(1) = rho*A*(x3d(1)-s1*(x2(1)-x2d(1))-xi*sign(s)-k*s-Omega(1))-f_c*sign(s)-rhoA*g(1);
        end

        f = sat(f);

    end
    function dy=glfdiff(y,t,gam)
        if strcmp(class(y),'function_handle')
            y=y(t);
        end
        h=t(2)-t(1); w=1; y=y(:); t=t(:);
        for j=2:length(t)
            w(j)=w(j-1)*(1-(gam+1)/(j-1));
        end
        for i=1:length(t)
            dy(i)=w(1:i)*[y(i:-1:1)]/h^gam;
        end
        dy = dy(end);
    end

end










