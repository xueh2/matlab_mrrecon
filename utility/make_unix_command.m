function gt_command = make_unix_command(user, gt_host, src_command)
% gt_command = make_unix_command(user, gt_host, src_command)

gt_command = ['ssh -o ConnectTimeout=10 -o TCPKeepAlive=yes -o ServerAliveInterval=15 -o ServerAliveCountMax=3 ' user '@' gt_host ' "' src_command '"'];