m = 10;
xidxs = xvalidationIdx(10 * m, 10);
for k = 1:length(xidxs)
    trainIdx = xidxs{k,1};
    testIdx = xidxs{k,2};
    assert(length(trainIdx) == 9 * m, 'training set size');
    assert(length(testIdx) == m, 'testing set size');
    disp(testIdx)
end

xidxs = xvalidationIdx(5, 2, false);
assert(all(xidxs{1,2} == [1 2 3]), 'I know the answer');
assert(all(xidxs{2,2} == [3 4 5]), 'I know the answer');

xidxs = xvalidationIdx(5, 3, false);
assert(all(xidxs{1,2} == [1 2]), 'I know the answer');
assert(all(xidxs{2,2} == [3 4]), 'I know the answer');
assert(all(xidxs{3,2} == [4 5]), 'I know the answer');

xidxs = xvalidationLPO((4:8)', (10:11)');
for k = 1:size(xidxs,1)
    disp(xidxs{k,2}');
    disp(xidxs{k,1}');
end
