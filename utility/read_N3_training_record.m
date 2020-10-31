function [inputs, targets, pred] = read_N3_training_record(train_dir, epoch, batch, plot_flag)
% read in N3 training record
% [inputs, targets, pred] = read_N3_training_record(train_dir, epoch, batch, plot_flag)

a = readNPY(fullfile(train_dir, ['inputs_epoch_' num2str(epoch) '__batch_' num2str(batch)]));
b = readNPY(fullfile(train_dir, ['targets_epoch_' num2str(epoch) '__batch_' num2str(batch)]));
c = readNPY(fullfile(train_dir, ['pred_epoch_' num2str(epoch) '__batch_' num2str(batch)]));

inputs = permute(a, [3 4 1 2]);
targets = permute(b, [3 4 1 2]);
pred = permute(c, [3 4 1 2]);

if(plot_flag)
    figure; imagescn(inputs, [], [], 12);
    figure; imagescn(targets, [], [], 12);
    figure; imagescn(pred, [], [], 12);
end