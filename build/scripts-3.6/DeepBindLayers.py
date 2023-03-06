from keras import backend as K
from keras.layers import Conv1D
from keras.utils.np_utils import conv_output_length


class RevCompConv1D(Conv1D):
    '''Like Convolution1D, except the reverse-complement filters with tied
    weights are added in the channel dimension. The reverse complement
    of the channel at index i is at index -i.

    # Example

    ```python
        # apply a reverse-complemented convolution 1d of length 20
        # to a sequence with 100bp input, with 2*64 output filters
        model = Sequential()
        model.add(RevCompConv1D(nb_filter=64, filter_length=20,
                                border_mode='same', input_shape=(100, 4)))
        # now model.output_shape == (None, 100, 128)

        # add a new reverse-complemented conv1d on top
        model.add(RevCompConv1D(nb_filter=32, filter_length=10,
                                border_mode='same'))
        # now model.output_shape == (None, 10, 64)
    ```

    # Arguments
        nb_filter: Number of non-reverse complemented convolution kernels
            to use (half the dimensionality of the output).
        filter_length: The extension (spatial or temporal) of each filter.
        init: name of initialization function for the weights of the layer
            (see [initializations](../initializations.md)),
            or alternatively, Theano function to use for weights initialization.
            This parameter is only relevant if you don't pass a `weights` argument.
        activation: name of activation function to use
            (see [activations](../activations.md)),
            or alternatively, elementwise Theano function.
            If you don't specify anything, no activation is applied
            (ie. "linear" activation: a(x) = x).
        weights: list of numpy arrays to set as initial weights
            (reverse-complemented portion should not be included as 
            it's applied during compilation)
        border_mode: 'valid', 'same' or 'full'. ('full' requires the Theano backend.)
        subsample_length: factor by which to subsample output.
        W_regularizer: instance of [WeightRegularizer](../regularizers.md)
            (eg. L1 or L2 regularization), applied to the main weights matrix.
        b_regularizer: instance of [WeightRegularizer](../regularizers.md),
            applied to the bias.
        activity_regularizer: instance of [ActivityRegularizer](../regularizers.md),
            applied to the network output.
        W_constraint: instance of the [constraints](../constraints.md) module
            (eg. maxnorm, nonneg), applied to the main weights matrix.
        b_constraint: instance of the [constraints](../constraints.md) module,
            applied to the bias.
        bias: whether to include a bias
            (i.e. make the layer affine rather than linear).
        input_dim: Number of channels/dimensions in the input.
            Either this argument or the keyword argument `input_shape`must be
            provided when using this layer as the first layer in a model.
        input_length: Length of input sequences, when it is constant.
            This argument is required if you are going to connect
            `Flatten` then `Dense` layers upstream
            (without it, the shape of the dense outputs cannot be computed).

    # Input shape
        3D tensor with shape: `(samples, steps, input_dim)`.

    # Output shape
        3D tensor with shape: `(samples, new_steps, nb_filter)`.
        `steps` value might have changed due to padding.
    '''

    def get_output_shape_for(self, input_shape):
        length = conv_output_length(input_shape[1],
                                    self.filter_length,
                                    self.border_mode,
                                    self.subsample[0])
        return (input_shape[0], length, 2*self.nb_filter)

    def call(self, x, mask=None):
        #create a rev-comped W. The last axis is the output channel axis.
        #dim 1 is dummy axis of size 1 (see 'build' method in Convolution1D)
        #Rev comp is along both the length (dim 0) and input channel (dim 2)
        #axes; that is the reason for ::-1, ::-1 in the first and third dims.
        #The rev-comp of channel at index i should be at index -i
        #This is the reason for the ::-1 in the last dim.
        rev_comp_W = K.concatenate([self.W, self.W[::-1,:,::-1,::-1]],axis=-1)
        if (self.bias):
            rev_comp_b = K.concatenate([self.b, self.b[::-1]], axis=-1)
        x = K.expand_dims(x, 2)  # add a dummy dimension
        output = K.conv2d(x, rev_comp_W, strides=self.subsample,
                          border_mode=self.border_mode,
                          dim_ordering='tf')
        output = K.squeeze(output, 2)  # remove the dummy dimension
        if self.bias:
            output += K.reshape(rev_comp_b, (1, 1, 2*self.nb_filter))
        output = self.activation(output)
        return output