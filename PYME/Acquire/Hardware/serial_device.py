import serial
import time
import threading
import logging
logger = logging.getLogger(__name__)

try:
    import Queue
except ImportError:
    import queue as Queue


class SerialDevice(object):
    def __init__(self, com_port, name, timeout=1):
        self.name = name
        self.timeout = timeout
        # initialize and configure the serial port without opening it
        self.com_port = serial.Serial(timeout=timeout)
        self.com_port.port = com_port
        self.lock = threading.Lock()
        self.is_on = False
        self.turn_on()

    def turn_on(self):
        self.is_on = True
        self.com_port.open()

    def query(self, command, lines_expected=1):
        with self.lock:
            self.com_port.reset_input_buffer()
            self.com_port.write(command)
            reply = [self.com_port.readline() for line in range(lines_expected)]
            self.com_port.reset_input_buffer()
        return reply

    def close(self):
        logger.debug('Shutting down %s' % self.name)
        # stop polling
        self.is_on = False
        self.com_port.close()



# class SerialDevice(object):
#     def __init__(self, com_port, name, timeout=1):
#         self.name = name
#         self.com_port = com_port
#         self.timeout = timeout
#
#         # initialize and configure the serial port without opening it
#         self.port = serial.Serial(timeout=timeout)
#         self.port.port = self.com_port
#
#         self.send_queue = Queue.Queue()
#         self.reply_queue = Queue.Queue()
#
#         self.is_on = False
#         # self.polling_thread = None
#         # self.lock = threading.Lock()
#         # self.turn_on()
#
#     def _purge_replies(self):
#         try:
#             while True:
#                 self.reply_queue.get_nowait()
#         except:
#             pass
#
#     def _command(self, command):
#         self.send_queue.put(command)
#
#     def query(self, command, replies_expected=1):
#         with self.lock:  # make sure no one takes our reply or sends commands before we're finished
#             self._purge_replies()
#             self.port.write(command)
#             reply = []
#             for ind in range(replies_expected):
#                 reply.append(self.reply_queue.get(timeout=self.timeout))
#
#             # if flush_after:
#             #     with self._internal_lock:
#             #         logger.debug('HERE')
#             #         with serial.Serial(self.com_port, timeout=self.timeout) as port:
#             #             port.flush()
#             self._purge_replies()
#             if replies_expected == 1:
#                 return reply[0]
#             return reply
#
#     def _poll(self):
#         while self.is_on:
#             logger.debug('here')
#             try:
#                 command = self.send_queue.get(False)
#                 logger.debug('command %s' % command)
#                 self.port.write(command)
#
#             except Queue.Empty:
#                 pass
#
#             listen = True
#             while listen == True:
#                 reply = self.port.readline()
#                 if reply != b'':
#                     self.reply_queue.put(reply)
#                 else:
#                     listen = False
#             time.sleep(.05)
#
#     def turn_on(self):
#         self.is_on = True
#
#         self.port.open()
#         # make sure we are polling
#         # if not self.polling_thread or not self.polling_thread.is_alive():
#         #     self.polling_thread = threading.Thread(target=self._poll)
#         #     self.polling_thread.start()
#
#     def close(self):
#         logger.debug('Shutting down %s' % self.name)
#         # stop polling
#         self.is_on = False
#         try:
#             self.polling_thread.join()
#         except:
#             pass
#
# class ContextConnectionSerialDevice(SerialDevice):
#     def _poll(self):
#         while self.is_on:
#             with self._internal_lock:
#                 with serial.Serial(self.com_port, timeout=self.timeout) as port:
#                     logger.debug('here')
#                     try:
#                         command = self.send_queue.get(False)
#                         logger.debug('command %s' % command)
#                         port.write(command)
#
#                     except Queue.Empty:
#                         pass
#
#                     listen = True
#                     while listen == True:
#                         reply = port.readline()
#                         if reply != b'':
#                             self.reply_queue.put(reply)
#                         else: listen = False
#             time.sleep(.05)
#
#     def turn_on(self):
#         self.is_on = True
#         # make sure we are polling
#         if not self.polling_thread or not self.polling_thread.is_alive():
#             self.polling_thread = threading.Thread(target=self._poll)
#             self.polling_thread.start()