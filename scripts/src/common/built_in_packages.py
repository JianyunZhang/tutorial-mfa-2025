import itertools as it
import warnings
from collections import defaultdict, abc
import gzip
import pickle
import pathlib
import os
import copy
import multiprocessing as mp
import argparse
import enum

# ValueEnum 是一个继承自 enum.Enum 的自定义枚举类，提供了更方便的字符串处理功能
class ValueEnum(enum.Enum):
    # 重写了 __str__ 方法
    def __str__(self):
        # 直接返回枚举值 (self.value)，而不是枚举成员的名称
        # 使得在打印或字符串转换时更直观
        return self.value

    # 添加了 startswith 方法
    # 检查枚举值是否以指定的子字符串开头
    # 直接代理到枚举值的 startswith 方法
    def startswith(self, substr):
        return self.value.startswith(substr)
