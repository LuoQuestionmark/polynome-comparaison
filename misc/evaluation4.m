function ret = evaluation4(X, tau)
  % arg1 = X(1);
  % arg2 = X(2);
  % [ A, b ] = alpha_generator2(arg1, arg2, tau);

  % Erreur = A * X(3 : 13) - b

  % % % % % % % % % % % % % % % % % % % % % % % % % %
  % % % % % % % % % % % % %
  % c_min = 1;
  % c_max = 1.5;
  a_min = 1;
  den2 = 2;
  b_min = -1;
  den0 = 1;
  a_max = 2;
  b_max = 1;
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  c_min = 1.25 - 0.25 * tau;
  c_max = 1.25 + 0.25 * tau;
  % c_min = 1;
  % den2 = 2;
  den0 = 1;
  % b_min = -tau;
  b_max = tau;
  % a_min = 1.25 - 0.25 * tau;
  a_max = 1.25 + 0.25 * tau;

  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
  % Num et denum du 1 système
  N1 = [c_min 0 a_min];
  D1 = [den2 b_min den0];

  % % Num et denum du 2 système
  N2 = [c_min 0 a_max];
  D2 = [den2 b_min den0];

  % % Num et denum du 3 système
  N3 = [c_min 0 a_min];
  D3 = [den2 b_max den0];

  % % Num et denum du 4 système

  N4 = [c_min 0 a_max];
  D4 = [den2 b_max den0];

  Nc = [0 X(11) 0];
  Dc = [X(12) 0 X(13)];

  % Stabilité des systèmes

  P1 = conv(N1, Nc) + conv(D1, Dc);
  RacP1 = roots(P1)

  P2 = conv(N2, Nc) + conv(D2, Dc);
  RacP2 = roots(P2)

  P3 = conv(N3, Nc) + conv(D3, Dc);
  RacP3 = roots(P3)

  P4 = conv(N4, Nc) + conv(D4, Dc);
  RacP4 = roots(P4)

  if any (real(RacP1) >= 0) || any(real(RacP2) >= 0) || any(real(RacP3) >= 0) || any(real(RacP4) >= 0) ret = false;
  else ret = true;
  end

end

