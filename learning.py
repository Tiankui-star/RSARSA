# from sympy import isprime
# #
# def uint32_array_to_int(arr):
#     """
#     将 uint32_t 数组（低位在前，高位在后）转换为大整数
#     """
#     n = 0
#     for i in reversed(range(len(arr))):  # 从高位开始
#         n = (n << 32) + arr[i]
#     return n
#
# def check_prime(uint32_arr):
#     n = uint32_array_to_int(uint32_arr)
#     print("Number:", n)
#     prime_flag = isprime(n)
#     print("Is prime?", prime_flag)
#     return n
# # # 示例：两个32位数组表示64位整数
# prime1 = [757897211,4128549411,3339880654,130291540,756149348,1356868423,2413761437,3391807815,3679400074,484489568,2746056864,2667321572,1959813599,1976730658,1777281610,4168282924]  # 低位在前，高位在后
# prime2 = [478585197,103485950,2904745936,2847268253,296791683,1193118780,1892927142,2044557110,855224915,2341389853,1166781217,3371591298,1400981518,532120207,1563684065,3426311462]  # 低位在前，高位在后
# t1=uint32_array_to_int(prime1)
# t2=uint32_array_to_int(prime2)

# from sympy import randprime
#
# # 生成64位素数
# p = randprime(2**63, 2**64)
# print(p)
# import struct

# def toarray(a):
#     arr=[]
#     while a>=0:
#         arr.append(a)
# def is_modular_inverse(a, b, modulus):
#     product = (a * b) % modulus
#

#
#     # 检查结果是否为1
#     return product == 1
#
#
# def int_to_uint32_array_le(num):
#     """
#     将大整数转换为小端序uint32数组
#     """
#     if num == 0:
#         return [0]
#
#     # 计算需要的字节数
#     byte_length = (num.bit_length() + 7) // 8
#     if byte_length == 0:
#         return [0]
#
#     # 转换为字节序列（小端序）
#     bytes_data = num.to_bytes(byte_length, byteorder='little', signed=False)
#
#     # 在后面填充0字节，使总字节数为4的倍数
#     padding_length = (4 - (len(bytes_data) % 4)) % 4
#     bytes_data = bytes_data + b'\x00' * padding_length
#
#     # 将字节序列转换为uint32数组（小端序）
#     uint32_array = []
#     for i in range(0, len(bytes_data), 4):
#         chunk = bytes_data[i:i + 4]
#         # 使用小端序解析
#         value = struct.unpack('<I', chunk)[0]
#         uint32_array.append(value)
#
#     return uint32_array
# mod=[4282406776,3798388802,3832662419,1331614735,3774489143,1646500511,2216546241,320687689,1463002823,4279859043,945193140,1632572325,2161928674,966202914,3434735367,2137567521,2989042459,2331711727,2955566397,2694762471,3692622802,1316537001,2811698907,2065487660,1682367034,2030934935,3805817619,4229741118,2519123284,1688665613,2267647529,3325248966]
# bb=[1002453145,2710426125,3490538576,2630149009,2947606774,3339188899,3708042018,3491890554,2209874560,3021965692,106305291,803725193,3925346375,3909576876,2475120960,119906250,670860622,750819242,2417715586,1491024699,2136699031,599231748,2071979722,325911475,1088590434,2577360045,3220523276,715730231,1882666084,1850601390,962011072,2151631684]
# mu=[1223921887,3735456868,1487354418,14207234,532462879,4196487715,2228267524,1462085319,1702660517,2810771169,563063926,3376517900,1227756496,3475053780,2480733746,1142227316,2989042461,2331711727,2955566397,2694762471,3692622802,1316537001,2811698907,2065487660,1682367034,2030934935,3805817619,4229741118,2519123284,1688665613,2267647529,3325248966]
# mumu=uint32_array_to_int(mu)
# modulus = uint32_array_to_int(mod)
# a = 17
# b = uint32_array_to_int(bb)  # 5*7=35 ≡ 1 mod 17
# tty=fast_power_mod(131072,b,mumu)
# array=int_to_uint32_array_le(tty)
# print(array)
# if __name__ == "__main__":
#     # 参数：vector as comma-separated string, d, operation
#     nums = list(map(int, sys.argv[1].split(',')))
#     d = int(sys.argv[2])
#     operation = sys.argv[3]
#
#     if operation == "add":
#         print(add(nums, d))
#     elif operation == "multiply":
#         print(multiply(nums, d))
#     elif operation == "max_plus":
#         print(max_plus_d(nums, d))
#     else:
#         print("Invalid operation")
import sys
def fast_power_mod(base, exponent, modulus):
    """
    快速幂算法计算 base^exponent mod modulus
    """
    result = 1
    base = base % modulus

    while exponent > 0:
        # 如果指数是奇数，乘上当前的base
        if exponent % 2 == 1:
            result = (result * base) % modulus

        # 指数右移一位，base平方
        exponent = exponent // 2
        base = (base * base) % modulus

    return result
def uint32_array_to_int(arr):
    """
    将 uint32_t 数组（低位在前，高位在后）转换为大整数
    """
    n = 0
    for i in reversed(range(len(arr))):  # 从高位开始
        n = (n << 32) + arr[i]
    return n
def int_to_uint32_array_le(n):
    """
    将大整数 n 转为小端 uint32 数组
    低位在前，高位在后
    """
    arr = []
    while n > 0:
        arr.append(n & 0xFFFFFFFF)  # 取低 32 位
        n >>= 32                     # 右移 32 位
    return arr
def jiemi(list1, list2, list3):
    # 把字符串转换成列表
    nums1 = list(map(int, list1.split(',')))
    nums2 = list(map(int, list2.split(',')))
    nums3 = list(map(int, list3.split(',')))

    # 示例处理：返回每个列表的和
    base = uint32_array_to_int(nums1)
    d = uint32_array_to_int(nums2)
    mod = uint32_array_to_int(nums3)

    rsu=fast_power_mod(base, d, mod)
    print (int_to_uint32_array_le(rsu))


if __name__ == "__main__":
    list1 = sys.argv[1]
    list2 = sys.argv[2]
    list3 = sys.argv[3]
    operation = sys.argv[4]

    if operation == "jiemi":
        jiemi(list1, list2, list3)
    else:
        print("Unknown operation")