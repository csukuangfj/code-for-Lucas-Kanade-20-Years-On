% lk20style.m
% Line styles for all LK20 figures

% Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
% $Id: lk20style.m,v 1.1.1.1 2003/08/20 03:07:36 iainm Exp $

% Figure text
font_size = 16;

% Forwards additive
ls.fa = 'r-x';
lw.fa = 5;
ms.fa = 14;

% Forwards compositional
ls.fc = 'g--s';
lw.fc = 3;
ms.fc = 12;

% Inverse additive
ls.ia = 'm-v';
lw.ia = 5;
ms.ia = 12;

% Inverse compositional
ls.ic = 'b--x';
lw.ic = 3;
ms.ic = 14;

% Inverse compositional - diagonal Hessian
ls.ic_d = 'b--*';
lw.ic_d = 3;
ms.ic_d = 12;

% Inverse compositional - steepest descent
ls.ic_sd = 'c-^';
lw.ic_sd = 3;
ms.ic_sd = 12;

% Inverse compositional - Newton
ls.ic_nt = 'g-p';
lw.ic_nt = 3;
ms.ic_nt = 12;

% Inverse compositional - Newton, diagonal Hessian
ls.ic_nt_d = 'g-*';
lw.ic_nt_d = 3;
ms.ic_nt_d = 12;

% Inverse compositional - Levenburg-Marquardt
ls.ic_lm = 'r-+';
lw.ic_lm = 2;
ms.ic_lm = 12;

% Inverse compositional - diagonal Hessian 2
ls.ic_d2 = 'b->';
lw.ic_d2 = 3;
ms.ic_d2 = 12;

% Inverse compositional - Newton, diagonal Hessian 2
ls.ic_nt_d2 = 'g->';
lw.ic_nt_d2 = 3;
ms.ic_nt_d2 = 12;

% Inverse compositional - reparameterised
ls.ric = 'm--x';
lw.ric = 2;
ms.ric = 11;

% Inverse compositional - diagonal Hessian - reparameterised
ls.ric_d = 'm--*';
lw.ric_d = 2;
ms.ric_d = 11;

% Inverse compositional - steepest descent - reparameterised
ls.ric_sd = 'g-^';
lw.ric_sd = 2;
ms.ric_sd = 11;
