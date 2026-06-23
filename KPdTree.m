addpath('qclab')
addpath('util')

% parameters
d = 5;
n = 4;

IMAG = 1;

% generate random normalized psi
[psi, ~] = qr(randn(d^n) + 1i * IMAG * randn(d^n));
psi = psi(:, 1);

% generate the tree with psi on the leaves
% and the (phased) norms of the children in the
% internal nodes
t = cell(1, n+1);

for i = 1:n
  t{i} = cell(1, d^(n-i));
  for j = 1:d^(n-i)
    if i == 1
      start = (j-1)*d;
      end_ = start + d - 1;
      t{i}{j} = psi(start+1:end_+1);
    else
      for k = 1:d
        c = t{i-1}{(j-1)*d + k}; %children
        t{i}{j}(k) = norm(c) * c(1) / abs(c(1));
      end
    end
  end
end

c = t{n}{1};
t{n+1} = norm(c) * c(1) / abs(c(1));

% replace the elements with the vector describing the
% householder reflector
e1 = zeros(d,1);
e1(1) = 1;
for i = 1:n
  for j = 1:d^(n-i)
    o = reshape(t{i}{j}, [], 1);
    t{i}{j} = o - sqrt(o'*o) * o(1)/abs(o(1)) * e1;
    t{i}{j} = t{i}{j} / norm(t{i}{j});
  end
end

cir = qclab.QCircuit(n, 0, d);

% add initial phase gate with element from root
phase = t{n+1};
phaseGate = qclab.qgates.Phase(0, real(phase), imag(phase));
phaseGate = qclab.qgates.qudit.SubspaceGate(phaseGate, [1,0], 0);
cir.push_back(phaseGate);

% foreach |v> in the tree, generate the householder
% and apply it in the circuit, with its corresponding
% controls
for i = n:-1:1
  for j = 1:numel(t{i})
    v = t{i}{j};
    W = eye(d) - 2 * (v*v');
    target = n-i;
    QHR = qclab.qgates.MatrixGate(target, W);
    if i < n
      ctrl = 0:n-i-1;
      ctrlStates = dec2base_(j-1, d, n-i);
      QHR = qclab.qgates.MControlledGate(QHR, ctrl, target, ctrlStates);
    end
    cir.push_back(QHR);
  end
end

% test whether the state of the circuit corresponds
% to the original psi
norm(cir.simulate(repmat('0', 1, n)).states - psi)

