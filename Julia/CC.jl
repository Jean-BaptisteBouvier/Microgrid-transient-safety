using Convex, SCS, MAT, Dates, CSDP

@time begin
cd("/Users/nand311/Documents/MATLAB/unc_analysis")
# Loading the open loop system matrices
file_read = matopen("openloop_sys_matrices.mat")
Ao = read(file_read, "A_o")
Bo = read(file_read, "B")
ng = read(file_read, "ng")
close(file_read)

ng = convert(Int,ng);

row_Ao, column_Ao = size(Ao);
row_Bo, column_Bo = size(Bo);

Q = Semidefinite(row_Bo,row_Bo)
M = Variable(column_Bo,row_Bo)
Alphai = Variable(column_Bo-ng, Positive())
Betai = Variable(column_Bo-ng, Positive())

expr1 = 0;

for i = ng+1:column_Bo
   expr1 = expr1 + Bo[:,i]*Alphai[i-ng]*Bo[:,i]';
end

# objective = norm(sum(Betai)-sum(Alphai));
objective = sum(Betai);
# objective = -sum(Alphai);

problem = minimize(objective);
constraints = ((-(Ao*Q + Q*Ao' + expr1 - Bo*M - (Bo*M)') - 1e-8*eye(row_Bo)) in :SDP);
#
problem.constraints += constraints;

for i = ng+1:column_Bo
    problem.constraints += (([Betai[i-ng] M[i,:]; M[i,:]' Q] - 1*1e-8*eye(row_Bo+1)) in :SDP);
end
# problem.constraints += (trace(Q) == 1);
solve!(problem, CSDPSolver())
# Alphai.value
# Betai.value

K = M.value*inv(Q.value);
Sigma = Alphai.value./Betai.value;

CurDir = pwd();
FolderName = Dates.format(now(),"YYYYudd_HH-MM");
ResDir = @sprintf("Dat_%s",FolderName);
mkdir(ResDir)

file_write = matopen(@sprintf("%s/%s/scenario.mat",CurDir,ResDir), "w")
write(file_write, "alpha", Alphai.value)
write(file_write, "beta", Betai.value)
write(file_write, "Q", Q.value)
write(file_write, "M", M.value)
write(file_write, "K", K)
write(file_write, "Sigma", Sigma)
close(file_write)

end
