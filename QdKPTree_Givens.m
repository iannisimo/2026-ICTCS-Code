addpath('qclab')
addpath('util')

% parameters
d = 3;
n = 3;


% generate random normalized psi
[psi, ~] = qr(randn(d^n));
psi = psi(:, 1);

% generate the tree with psi on the leaves
% and the (phased) norms of the children in the
% internal nodes
t = cell(1, n+1);
s = cell(1, 1);

for i = 1:n
  t{i} = cell(1, d^(n-i));
  if (i == 1)
    s{i} = cell(1, d^(n-i));
  end
  for j = 1:d^(n-i)
    if i == 1
      start = (j-1)*d;
      end_ = start + d - 1;
      t{i}{j} = power(psi(start+1:end_+1),2);
      s{i}{j} = sign(real(psi(start+1:end_+1))) + 1i * sign(imag(psi(start+1:end_+1)));
    else
      for k = 1:d
        t{i}{j}(k) = sum(t{i-1}{(j-1)*2 + k});
      end
    end
  end
end

t{n+1} = sqrt(sum(power(psi, 2))) / norm(psi);

% convert into binary tree for givens
for i = 1:n
  for j = 1:d^(n-i)
    dary = t{i}{j};
    dary = reshape(dary, 1, []);
    dary = [dary, zeros(1, (power(2, ceil(log2(d))) - numel(dary)))];
    sub_depth = ceil(log2(d));
    t{i}{j} = cell(1, sub_depth);
    for k = 1:sub_depth
      D = numel(dary);
      elems = D / (2^k);
      t{i}{j}{k} = cell(1, elems);
      for l = 1:elems
        if k == 1
          t{i}{j}{k}{l} = dary((l-1)*2 + 1:l*2);
        else
          for m = 1:2
            t{i}{j}{k}{l}(m) = sum(t{i}{j}{k-1}{(l-1) * 2 + m});
          end
        end
      end
    end
  end
end

cir = qclab.QCircuit(n, 0, d);


for i = n:-1:1
  for j = 1:numel(t{i})
    % set controls and target
    if i < n
      ctrl = 0:n-i-1;
      ctrlStates = dec2base_(j-1, d, n-i);
    end
    target = n-i;
    sub_depth = numel(t{i}{j});
    for k = sub_depth:-1:1
      for l = 1:numel(t{i}{j}{k})
        old = 1;
        if (k < sub_depth)
          old = t{i}{j}{k+1}{ceil(l/2)}(mod(l-1, 2) + 1);
        elseif i < n
          child_idx = mod(j-1, d) + 1;
          old = t{i+1}{ceil(j/d)}{1}{ceil(child_idx/2)}(mod(child_idx-1, 2) + 1);
        end
        phi = t{i}{j}{k}{l};
        if(any(phi == 0)), continue; end
        ry = qclab.qgates.RotationY(target, sqrt(phi(1)) / sqrt(old), sqrt(phi(2)) / sqrt(old), true);
        % givens subspace calculation
        leaf_l = (l-1) * 2;
        leaf_r = leaf_l + 1;
        sub_l = leaf_l * 2^(k-1);
        sub_r = leaf_r * 2^(k-1);
        ry = qclab.qgates.qudit.SubspaceGate(ry, [sub_l, sub_r], target);
        if i < n
          ry = qclab.qgates.MControlledGate(ry, ctrl, target, ctrlStates);
        end
        cir.push_back(ry);
      end
    end
  end
end

% In this demo, signs are fixed via additional phase gates.
% This operation should be integrated in the rotations as
% for the qubit algorithm, however this method functionally
% has the same result.
for i = 1:numel(s{1})
  for j = 1:numel(s{1}{i})
    if sign(s{1}{i}(j)) == -1
      p = qclab.qgates.Phase(n-1, pi);
      p = qclab.qgates.qudit.SubspaceGate(p, [mod(j, d), j-1], n-1);
      if n > 0
        ctrl = 0:n-2;
        ctrlStates = dec2base_(i-1, d, n-1);
        p = qclab.qgates.MControlledGate(p, ctrl, n-1, ctrlStates);
      end
      cir.push_back(p);
    end
  end
end

% test whether the state of the circuit corresponds
% to the original psi
normdiff = norm(cir.simulate(repmat('0', 1, n)).states - psi);
fprintf('Norm of the difference between psi and the simulated state: %.2f\n', normdiff);
