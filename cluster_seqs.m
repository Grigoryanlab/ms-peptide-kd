% read sequences and the substitution matrix
F = fastaread('peptides.fasta');
alpha = 'ACDEFGHIKLMNPQRSTVWY';
B = blosum(65, 'Order', alpha);

% encode sequences via substitution-matrix rows for each amino acid
S = zeros(length(F), length(F(1).Sequence)*length(alpha));
G = cell(length(F), 1); % groupings by the P(0) amino acid
for i = 1:length(F)
    for j = 1:length(F(i).Sequence)
        aaj = find(alpha == F(i).Sequence(j), 1);
        S(i, length(alpha)*(j-1)+1:length(alpha)*j) = B(aaj, :);
    end
    G{i} = F(i).Sequence(end);
end

% do a t-SNE plot
Y = tsne(S);
gscatter(Y(:,1), Y(:,2), G);

% make it pretty
h = gscatter(Y(:,1), Y(:,2), G, [], '.', 12);
set(gca, 'FontSize', 16);
axis square;
xlabel('t-SNE dimension 1');
ylabel('t-SNE dimension 2');
set(gca, 'XTick', [], 'YTick', []);
print(gcf, '-dpng', '-r300', 'tsne.png');
