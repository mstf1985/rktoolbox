
fprintf('\n Please wait while running all tests:\n\n')

tests = dir('test_*.m');
failed = 0;
for j = 1:length(tests),

    fn = tests(j).name; 
    fn = fn(1:end-2);
    fprintf('   %4s %18s ... ',['[' num2str(j) ']'],fn)
    try
        check = feval(fn);
    catch e,
        check = 0;
    end
    if ~all(check),
        fprintf('failed.\n')
        failed = failed + 1;
    else
        fprintf('passed.\n')
    end

end

fprintf('   -----------------------------------\n')
if failed,
    fprintf('   %2d tests failed. \n\n',failed)
else
    fprintf('   All tests passed. \n\n')
end