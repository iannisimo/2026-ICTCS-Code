addpath('qclab')
addpath('util')

% parameters
n = 5;


% generate random normalized psi
[psi, ~] = qr(randn(2^n));
psi = psi(:, 1);

% generate the tree with the square of the elements
% of psi on the leaves, and the sum of the children
% on the internal nodes; keep sign in a separate ordered
% list corresponding to the leaves of t
t = cell(1, n+1);
s = cell(1,1);

for i = 1:n
  t{i} = cell(1, 2^(n-i));
  if(i == 1)
    s{i} = cell(1, 2^(n-i));
  end
  for j = 1:2^(n-i)
    if i == 1
      start = (j-1)*2;
      end_ = start + 2 - 1;
      t{i}{j} = power(psi(start+1:end_+1),2);
      s{i}{j} = sign(real(psi(start+1:end_+1))) + 1i * sign(imag(psi(start+1:end_+1)));
    else
      for k = 1:2
        t{i}{j}(k) = sum(t{i-1}{(j-1)*2 + k});
      end
    end
  end
end

t{n+1} = sqrt(sum(power(psi, 2))) / norm(psi);

cir = qclab.QCircuit(n, 0, 2);

% foreach pair in the tree, generate the
% corresponding ry. For the last step,
% add the sign data
for i = n:-1:1
  for j = 1:numel(t{i})
    old = 1;
    if i < n
      old = t{i+1}{ceil(j/2)}(mod(j-1, 2) + 1);
    end
    phi = t{i}{j};
    cs = 1;
    ss = 1;
    if i == 1
      cs = s{i}{j}(1);
      ss = s{i}{j}(2);
    end
    ry = qclab.qgates.RotationY(n-i, cs*sqrt(phi(1))/sqrt(old), ss*sqrt(phi(2))/sqrt(old), true);
    if i < n
      ctrl = 0:n-i-1;
      target = n-i;
      ctrlStates = dec2base_(j-1, 2, n-i);
      ry = qclab.qgates.MControlledGate(ry, ctrl, target, ctrlStates);
    end
    cir.push_back(ry);
  end
end

% test whether the state of the circuit corresponds
% to the original psi
normdiff = norm(cir.simulate(repmat('0', 1, n)).states - psi);
fprintf('Norm of the difference between psi and the simulated state: %.2f\n', normdiff);
