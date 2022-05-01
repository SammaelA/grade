import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = "2"
import tensorflow as tf
import pandas as pd
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Input, Flatten, Dropout, Activation, BatchNormalization, Dense
from keras_preprocessing.image import ImageDataGenerator
from tensorflow import keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.applications import ResNet50
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
  tf.config.experimental.set_memory_growth(gpu, True)
dir = "../saves/NE_dataset/"
sess = tf.compat.v1.InteractiveSession()

df = pd.read_csv(dir +"dataset.csv")
y_columns = []
features_count = 1 #df.shape[1] - 1
for i in range(features_count):
    y_columns.append("par_"+str(i))
y_columns = ['par_2']
print(df.shape[1])
print(y_columns)
datagen=ImageDataGenerator(rescale=1./255.,validation_split=0.2)
#print(df)
train_df = df[:int(0.9*len(df))]
test_df = df[int(0.9*len(df)):]
B_SIZE = 8
train_generator=datagen.flow_from_dataframe(dataframe=train_df,x_col="image_dir",y_col=y_columns,subset="training",batch_size=B_SIZE,
                                            shuffle=True, class_mode="raw", target_size=(256,256))

valid_generator=datagen.flow_from_dataframe(dataframe=train_df,x_col="image_dir",y_col=y_columns,subset="validation",batch_size=B_SIZE,
                                            shuffle=True, class_mode="raw", target_size=(256,256))

test_datagen=ImageDataGenerator(rescale=1./255.)
test_generator=test_datagen.flow_from_dataframe(dataframe=test_df,x_col="image_dir",y_col=None,batch_size=B_SIZE,
                                                shuffle=False, class_mode=None, target_size=(256,256))

base_model = ResNet50(include_top=False,
                      input_shape=(256, 256, 3),
                      weights=None,
                      pooling='avg')

model = Sequential()

#model.add(layers.AveragePooling2D(pool_size=(2, 2), input_shape=(256, 256, 3)))
model.add(base_model)
model.add(Dense(features_count, activation="sigmoid"))
model.compile(loss='mse', optimizer=keras.optimizers.Adam(learning_rate=0.001),metrics=['mae'])
#base_model.summary()
model.summary()

STEP_SIZE_TRAIN=train_generator.n//train_generator.batch_size
STEP_SIZE_VALID=valid_generator.n//valid_generator.batch_size
model.fit(train_generator,
          batch_size=B_SIZE,
          epochs=1,
          validation_data=valid_generator,
          steps_per_epoch=STEP_SIZE_TRAIN,
          validation_steps=STEP_SIZE_VALID)
r = model.predict(test_generator)
print(r)