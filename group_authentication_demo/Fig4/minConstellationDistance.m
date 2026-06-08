function [dmin, idxPair, pPair] = minConstellationDistance(S)
% 输入：
%   S      : 列向量/任意形状的复数星座点集合（将会被拉直）
% 输出：
%   dmin   : 最小星座距离（欧氏距离）
%   idxPair: 产生最小距离的两点在原始 S 中的下标 [i j]
%   pPair  : 这两点的取值 [S(i) S(j)]

    % 拉直 & 过滤非法值
    S = S(:);
    finiteMask = isfinite(S);
    S = S(finiteMask);
    if numel(S) < 2
        dmin   = NaN;
        idxPair = [];
        pPair   = [];
        warning('有效星座点少于2个，无法计算最小距离。');
        return;
    end

    % 先检查是否有重复点（重复则最小距离为0，直接返回）
    [u, ~, ic] = unique(S, 'stable');           % u是去重后的点，ic是S到u的映射
    counts = accumarray(ic, 1);
    dupClass = find(counts > 1, 1, 'first');    % 找到第一个出现重复的类
    if ~isempty(dupClass)
        % 找到重复类在原始（过滤后）索引中的两个位置
        allIdx = find(ic == dupClass);
        i1 = allIdx(1); i2 = allIdx(2);
        % 映射回原始 S 输入中的下标（考虑过滤前）
        originalIdx = find(finiteMask);
        idxPair = originalIdx([i1 i2]);
        pPair   = S([i1 i2]);
        dmin    = 0;
        return;
    end

    % 否则计算任意两点的欧氏距离最小值
    n = numel(u);
    % 利用复数差的模长就是二维欧氏距离
    D = abs(u - u.');                 % n x n
    D(1:n+1:end) = Inf;               % 对角线置为无穷，避免选到自己
    [dmin, lin] = min(D(:));
    [i, j] = ind2sub([n n], lin);

    % 把去重后的下标映射回原始S的下标
    % unique(...,'stable') 返回的 u(i) 对应原始第一个出现的位置
    % 找到原始里与 u(i), u(j) 匹配的下标（第一个即可）
    originalIdx = find(finiteMask);
    ii = find(S == u(i), 1, 'first');
    jj = find(S == u(j), 1, 'first');
    idxPair = originalIdx([ii jj]);
    pPair   = [u(i) u(j)];
end
