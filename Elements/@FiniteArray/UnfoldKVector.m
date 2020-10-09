function kvec = UnfoldKVector(kvec, integrationpath)
    % Unfold the integration path. We assume that integrationpath(1) is the leftmost waypoint, and
    % integrationpath(end) is the rightmost waypoint.
    % The path is unfolded by first converting the imaginary tails into real parts. This is shown in
    % the illustration below, where * are waypoints, and the paths are shown using dots
    % Old Integration path (dotted)             New Integration path (dotted)
    %                |                                         |
    %                |*...........*                            |*............*
    %                :             :                           :              :
    % --------------:+--------------*-          --------------:+---------------*............
    %              * |              :           .............* |
    %              : |              :                          |
    %              : |              :                          |
    %              V |              V                          |
    % This conversion should (hopefully) retain relative distances
    
    % Determine the tails to unfold.
    leftTail = logical((real(kvec) < 0) & (imag(kvec) < imag(integrationpath(1))));
    rightTail = logical((real(kvec) > 0) & (imag(kvec) < imag(integrationpath(end))));
    
    % Convert the relative imaginary distance to a relative real distance
    kvec(leftTail) = integrationpath(1) + imag(kvec(leftTail) - integrationpath(1));
    kvec(rightTail) = integrationpath(end) - imag(kvec(rightTail) - integrationpath(end));
end
